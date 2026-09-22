#!/usr/bin/env python3
"""Resource-limited real-data comparison against the unmodified original R.

Run from the parent repository, after tests/prepare-carnivora.py:
  python3 tests/benchmark-carnivora.py build/carnivora-benchmark

Linux developer helper; not a dependency of the native CLI. All runs are serial.
The 3,000-gene C++ sample is an optional scaling check, not an R comparison.
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
parser.add_argument("--sizes", nargs="+", type=int, default=[100, 250, 1000])
parser.add_argument("--repeats", type=int, default=3)
parser.add_argument("--cpp-only", action="store_true")
parser.add_argument("--limit-gib", type=float, default=3)
parser.add_argument("--timeout", type=float, default=600)
args = parser.parse_args()
if args.repeats < 1 or not 0 < args.limit_gib <= 4 or args.timeout <= 0:
    parser.error("positive repeats/timeout and 0 < limit-gib <= 4 required")
repo = Path(__file__).resolve().parents[1]
output = args.output.resolve()
output.relative_to(repo)
output.mkdir(parents=True, exist_ok=False)
(output / "tmp").mkdir()
source = repo / "validation-data/carnivora-full"
binary = repo / "build/phylter"
reference = repo / ".audit/reference-lib"
limit = int(args.limit_gib * 1024**3)
reserve = 2 * 1024**3


def available():
    values = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    return int(values["MemAvailable"].split()[0]) * 1024


def restrict():
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


env = os.environ.copy()
env.update(TMPDIR=str(output / "tmp"), R_LIBS=f"{reference}:{repo / '.audit/library'}",
           TZ="Europe/Athens", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
metadata = {
    "platform": platform.platform(), "python": platform.python_version(),
    "cpu": next((s.split(":", 1)[1].strip() for s in
        Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")), "unknown"),
    "reference_commit": "4d74241169be3882da9d5f08da47bd40604764eb",
    "reference_library": str(reference),
    "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
    "data": json.loads((source / "provenance.json").read_text()),
    "threads": 1, "address_space_limit_bytes": limit, "memory_reserve_bytes": reserve,
    "initial_available_bytes": available(), "timeout_seconds": args.timeout,
    "command_line": vars(args) | {"output": str(output)},
}
(output / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")
rows, summary = [], []
for count in args.sizes:
    input_path = source / f"sample-{count}.nwk"
    if not input_path.is_file():
        raise ValueError(f"No prepared sample: {input_path}")
    folder = output / str(count)
    folder.mkdir()
    errors = []
    for repeat in range(args.repeats):
        prefixes = {}
        order = ["C++"] if args.cpp_only else (["R", "C++"] if repeat % 2 == 0 else ["C++", "R"])
        for implementation in order:
            free = available()
            if free < limit + reserve:
                raise RuntimeError(f"Not starting: available {free/1024**3:.2f} GiB < limit + 2 GiB reserve")
            prefix = folder / f"{implementation.replace('+', 'p')}-{repeat+1}"
            prefixes[implementation] = prefix
            command = (["Rscript", "--vanilla", str(repo / "tests/benchmark-worker.R"),
                        str(reference), str(input_path), str(prefix)] if implementation == "R" else
                       [str(binary), "run", "--trees", str(input_path), "--out", str(prefix), "--threads", "1"])
            print(f"START {count} genes / {implementation} / repeat {repeat+1}; available {free/1024**3:.2f} GiB", flush=True)
            with Path(str(prefix) + ".stdout").open("w") as stdout, Path(str(prefix) + ".stderr").open("w") as stderr:
                start = time.monotonic()
                child = subprocess.Popen(["/usr/bin/time", "-f", "%e\t%U\t%S\t%M", "-o",
                    str(prefix) + ".time.tsv", *command], env=env, stdout=stdout, stderr=stderr,
                    start_new_session=True, preexec_fn=restrict)
                try:
                    while child.poll() is None:
                        if time.monotonic() - start > args.timeout:
                            raise RuntimeError("Stopped at wall-time limit")
                        if available() < reserve:
                            raise RuntimeError("Stopped to preserve 2 GiB available memory")
                        time.sleep(0.2)
                    if child.returncode != 0:
                        raise RuntimeError(f"Benchmark exited {child.returncode}; inspect {prefix}.stderr")
                except BaseException:
                    try:
                        os.killpg(child.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    child.wait()
                    raise
            elapsed, user, system, rss = Path(str(prefix) + ".time.tsv").read_text().split()
            rows.append({"genes": count, "implementation": implementation, "repeat": repeat+1,
                         "elapsed_seconds": float(elapsed), "user_seconds": float(user),
                         "system_seconds": float(system), "peak_rss_mib": int(rss)/1024,
                         "r_analysis_seconds": float(Path(str(prefix) + ".analysis-seconds.txt").read_text())
                         if implementation == "R" else None})
            print(f"DONE  {count} / {implementation}: {elapsed}s, {int(rss)/1024:.2f} MiB peak RSS", flush=True)
            with (output / "runs.tsv").open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
        if not args.cpp_only:
            for field in ("outliers", "discarded"):
                assert table(Path(str(prefixes["R"]) + f".{field}.tsv")) == table(Path(str(prefixes["C++"]) + f".{field}.tsv")), f"{count}: {field} differ"
            rs, cs = (table(Path(str(prefixes[impl]) + ".scores.tsv")) for impl in ("R", "C++"))
            assert [x["state"] for x in rs] == [x["state"] for x in cs], "Accepted states differ"
            error = max(abs(float(r["quality"]) - float(c["quality"])) for r, c in zip(rs, cs))
            assert error < 1e-10, f"Quality differs by {error}"
            errors.append(error)
            print(f"VERIFIED ordered outliers/discards + all {len(rs)} accepted states; score error {error:.3g}", flush=True)
    item = {"genes": count, "taxa": 53, "repeats": args.repeats,
            "outliers": len(table(Path(str(prefixes["C++"]) + ".outliers.tsv"))),
            "discarded_pairs": len(table(Path(str(prefixes["C++"]) + ".discarded.tsv"))),
            "accepted_states": len(table(Path(str(prefixes["C++"]) + ".scores.tsv"))),
            "max_score_error": max(errors) if errors else None}
    for impl in order:
        sample = [r for r in rows if r["genes"] == count and r["implementation"] == impl]
        item[impl] = {"median_elapsed_seconds": statistics.median(r["elapsed_seconds"] for r in sample),
                      "elapsed_range_seconds": [min(r["elapsed_seconds"] for r in sample), max(r["elapsed_seconds"] for r in sample)],
                      "median_peak_rss_mib": statistics.median(r["peak_rss_mib"] for r in sample),
                      "max_peak_rss_mib": max(r["peak_rss_mib"] for r in sample)}
        if impl == "R":
            item[impl]["median_analysis_seconds"] = statistics.median(r["r_analysis_seconds"] for r in sample)
    summary.append(item)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
