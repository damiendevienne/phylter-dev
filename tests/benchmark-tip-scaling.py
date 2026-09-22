#!/usr/bin/env python3
"""Bounded, serial benchmark across real genes and synthetic taxa.

Python/R are validation dependencies only. Run after expand-carnivora-tips.py.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import signal
import statistics
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
parser.add_argument("--cpp-repeats", type=int, default=3)
parser.add_argument("--r-repeats", type=int, default=1)
parser.add_argument("--timeout", type=float, default=180)
parser.add_argument("--data", type=Path, help="Expanded dataset directory (default: validation-data/carnivora-tip-scaling)")
parser.add_argument("--initial-only",action="store_true",help="C++ initialization only; requires --r-repeats 0")
parser.add_argument("--reference-results",type=Path,help="Full-analysis result directory used to verify initial-only scores")
args = parser.parse_args()
if args.cpp_repeats <= 0 or args.r_repeats < 0 or args.timeout <= 0:
    parser.error("positive cpp-repeats/timeout and nonnegative r-repeats required")
if args.initial_only and args.r_repeats:
    parser.error("--initial-only requires --r-repeats 0; validate scores against the full-analysis reference runs")
if args.initial_only and not args.reference_results:
    parser.error("--initial-only requires --reference-results")
repo = Path(__file__).resolve().parents[2]
root = repo / "phylter2"
data = args.data.resolve() if args.data else root / "validation-data/carnivora-tip-scaling"
out = args.output.resolve()
out.relative_to(repo)
out.mkdir(parents=True, exist_ok=False)
(out / "tmp").mkdir()
binary = root / "build/phylter"
reference = repo / ".audit/reference-lib"
limit, reserve = 3*1024**3, 2*1024**3
env = os.environ.copy()
env.update(TMPDIR=str(out / "tmp"), R_LIBS=f"{reference}:{repo / '.audit/library'}",
           TZ="Europe/Athens", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")


def available():
    return int(next(line.split()[1] for line in Path("/proc/meminfo").read_text().splitlines()
                    if line.startswith("MemAvailable:"))) * 1024


def restrict():
    resource.setrlimit(resource.RLIMIT_AS, (limit,limit))
    resource.setrlimit(resource.RLIMIT_CORE, (0,0))


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


provenance = json.loads((data / "provenance.json").read_text())
metadata = {"data": provenance, "reference_commit": "4d74241169be3882da9d5f08da47bd40604764eb",
    "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
    "platform": platform.platform(), "python": platform.python_version(), "threads": 1,
    "limit_bytes": limit, "reserve_bytes": reserve, "timeout_seconds": args.timeout,
    "cpp_repeats": args.cpp_repeats, "r_repeats": args.r_repeats,"initial_only":args.initial_only,
    "cpu": next(line.split(":",1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name"))}
(out / "environment.json").write_text(json.dumps(metadata,indent=2)+"\n")
if args.initial_only:
    full_root = args.reference_results.resolve()
    previous = json.loads((full_root / "environment.json").read_text())
    assert previous["data"]==provenance, "Reference input provenance differs"
    assert previous["binary_sha256"]==metadata["binary_sha256"], "Reference C++ executable differs"
    metadata["reference_results"] = str(full_root)
    (out / "environment.json").write_text(json.dumps(metadata,indent=2)+"\n")
rows, summaries = [], []


def save():
    with (out / "runs.tsv").open("w") as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (out / "summary.json").write_text(json.dumps(summaries,indent=2)+"\n")


for case in sorted(provenance["datasets"],key=lambda x:x["genes"]*x["taxa"]**2):
    genes, taxa = case["genes"], case["taxa"]
    folder = out / f"g{genes}-t{taxa}"
    folder.mkdir()
    input_path = data / case["file"]
    assert hashlib.sha256(input_path.read_bytes()).hexdigest()==case["sha256"]
    good = {"R": [], "C++": []}
    stopped = set()
    for repeat in range(max(args.cpp_repeats,args.r_repeats)):
        order = ["R","C++"] if repeat%2==0 else ["C++","R"]
        for impl in order:
            if impl in stopped or repeat >= (args.cpp_repeats if impl=="C++" else args.r_repeats):
                continue
            prefix = folder / f"{impl.replace('+','p')}-{repeat+1}"
            row = {"genes":genes,"taxa":taxa,"implementation":impl,"repeat":repeat+1,
                   "status":"pending","elapsed_seconds":None,"user_seconds":None,
                   "system_seconds":None,"peak_rss_mib":None,"accepted_states":None,
                   "outliers":None,"r_analysis_seconds":None}
            # Conservative screening, not a prediction: R retains more matrices
            # and full RV temporaries. A hard limit independently bounds allocations.
            estimate = (64*1024**2+48*genes*taxa**2 if impl=="C++" else
                        256*1024**2+96*genes*taxa**2+48*genes**2)
            if estimate > limit:
                row["status"] = "skipped_memory_screen"
            elif available() < limit+reserve:
                row["status"] = "skipped_available_memory"
            else:
                command = (["Rscript","--vanilla",str(root / "tests/benchmark-worker.R"),str(reference),str(input_path),str(prefix)]
                           if impl=="R" else [str(binary),"run","--trees",str(input_path),"--out",str(prefix),"--threads","1"])
                if args.initial_only:
                    command.append("--initial-only")
                print(f"START {genes} genes / {taxa} taxa / {impl} / {repeat+1}",flush=True)
                with Path(str(prefix)+".stdout").open("w") as stdout, Path(str(prefix)+".stderr").open("w") as stderr:
                    child = subprocess.Popen(["/usr/bin/time","-f","%e\t%U\t%S\t%M","-o",str(prefix)+".time.tsv",*command],
                        env=env,stdout=stdout,stderr=stderr,start_new_session=True,preexec_fn=restrict)
                    start = time.monotonic()
                    try:
                        while child.poll() is None:
                            if time.monotonic()-start > args.timeout or available() < reserve:
                                row["status"] = "timeout" if time.monotonic()-start > args.timeout else "stopped_memory_reserve"
                                os.killpg(child.pid,signal.SIGKILL)
                                child.wait()
                                break
                            time.sleep(.2)
                        if row["status"]=="pending":
                            row["status"] = "ok" if child.returncode==0 else f"failed_exit_{child.returncode}"
                    except BaseException:
                        try:
                            os.killpg(child.pid,signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        child.wait()
                        raise
                if row["status"]=="ok":
                    elapsed,user,system,rss = Path(str(prefix)+".time.tsv").read_text().split()
                    row.update(elapsed_seconds=float(elapsed),user_seconds=float(user),system_seconds=float(system),
                               peak_rss_mib=int(rss)/1024,accepted_states=len(table(Path(str(prefix)+".scores.tsv"))),
                               outliers=len(table(Path(str(prefix)+".outliers.tsv"))))
                    if impl=="R":
                        row["r_analysis_seconds"] = float(Path(str(prefix)+".analysis-seconds.txt").read_text())
                    if good[impl]:
                        for field in ("scores","outliers","discarded"):
                            assert Path(str(prefix)+f".{field}.tsv").read_bytes()==Path(str(good[impl][0])+f".{field}.tsv").read_bytes(), "Repeat results differ"
                    good[impl].append(prefix)
            if row["status"]!="ok":
                stopped.add(impl)
            rows.append(row)
            print(f"DONE {genes}/{taxa}/{impl}: {row['status']}; {row['elapsed_seconds']} s; {row['peak_rss_mib']} MiB",flush=True)
            save()
    item = {"genes":genes,"taxa":taxa,"comparison":"not_available"}
    if good["R"] and good["C++"]:
        rp,cp = good["R"][0],good["C++"][0]
        for field in ("outliers","discarded"):
            assert table(Path(str(rp)+f".{field}.tsv"))==table(Path(str(cp)+f".{field}.tsv")), f"{folder}: {field} differ"
        rs,cs = (table(Path(str(p)+".scores.tsv")) for p in (rp,cp))
        assert [r["state"] for r in rs]==[c["state"] for c in cs],f"{folder}: states differ"
        error = max(abs(float(r["quality"])-float(c["quality"])) for r,c in zip(rs,cs))
        assert error < 1e-10,f"{folder}: quality error {error}"
        item.update(comparison="pass",max_score_error=error)
        print(f"VERIFIED {genes}/{taxa}: ordered outliers, discards, all states; max score error {error:.3g}",flush=True)
    if args.initial_only and good["C++"]:
        scores = table(Path(str(good["C++"][0])+".scores.tsv"))
        assert len(scores)==1 and scores[0]["state"]=="0", "Expected only the initial state"
        assert not table(Path(str(good["C++"][0])+".outliers.tsv")), "Unexpected initial outliers"
        errors, references = [], []
        for impl in ("Cpp","R"):
            path = full_root / folder.name / f"{impl}-1.scores.tsv"
            if path.exists():
                expected = table(path)[0]
                errors.append(abs(float(scores[0]["quality"])-float(expected["quality"])))
                references.append(impl)
        assert "Cpp" in references and max(errors)<1e-10, "Initial score differs from full-analysis reference"
        item.update(comparison="initial_matches_full_"+"_and_".join(references),max_score_error=max(errors))
        print(f"VERIFIED initial {genes}/{taxa} against {references}; error {max(errors):.3g}",flush=True)
    for impl in ("R","C++"):
        sample = [r for r in rows if r["genes"]==genes and r["taxa"]==taxa and r["implementation"]==impl]
        completed = [r for r in sample if r["status"]=="ok"]
        item[impl] = {"statuses":[r["status"] for r in sample],"repeats":len(completed)}
        if completed:
            for field in ("elapsed_seconds","peak_rss_mib","accepted_states","outliers"):
                item[impl][field] = statistics.median(r[field] for r in completed)
            item[impl]["elapsed_range_seconds"] = [min(r["elapsed_seconds"] for r in completed),max(r["elapsed_seconds"] for r in completed)]
    summaries.append(item)
    save()
print(json.dumps(summaries,indent=2))
