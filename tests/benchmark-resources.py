#!/usr/bin/env python3
"""Developer-only Linux benchmark; Python and R are not native runtime dependencies.

Run from the parent repository:
python3 tests/benchmark-resources.py build/resource-benchmark

All runs are sequential, in fresh processes, with a single thread.
GNU time measures the child process, excluding this Python orchestrator.
"""
import argparse
import csv
import hashlib
import json
import os
import platform
from pathlib import Path
import statistics
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output")
parser.add_argument("--multipliers", type=int, nargs="+", default=[1, 4])
parser.add_argument("--repeats", type=int)
parser.add_argument("--baseline", type=Path, help="Optional preserved previous C++ executable")
args = parser.parse_args()
if any(multiplier < 1 for multiplier in args.multipliers) or (args.repeats is not None and args.repeats < 1):
    parser.error("multipliers and repeats must be positive")
repo = Path(__file__).resolve().parents[1]
destination = Path(args.output).resolve()
destination.relative_to(repo)  # Keep every generated file inside this repository.
destination.mkdir(parents=True, exist_ok=False)
(destination / "tmp").mkdir()
binary = repo / "build/phylter"
baseline = args.baseline.resolve() if args.baseline else None
implementations = ["R", "C++"] + (["C++ dense"] if baseline else [])
reference = repo / ".audit/reference-lib"
environment = os.environ.copy()
environment.update({
    "TMPDIR": str(destination / "tmp"),
    "R_LIBS": str(reference) + os.pathsep + str(repo / ".audit/library"),
    "TZ": "Europe/Athens",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
})
source = (repo / "tests/fixtures/carnivora.nwk").read_text()
metadata = {
    "platform": platform.platform(),
    "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
    "fixture_sha256": hashlib.sha256(source.encode()).hexdigest(),
    "reference_commit": "4d74241169be3882da9d5f08da47bd40604764eb",
    "threads": 1,
}
if baseline:
    metadata["baseline_sha256"] = hashlib.sha256(baseline.read_bytes()).hexdigest()
(destination / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")
rows = []
report = []


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


for multiplier in args.multipliers:
    name = "carnivora" if multiplier == 1 else f"carnivora_{multiplier}x_genes"
    repeats = args.repeats or (5 if multiplier == 1 else 3)
    folder = destination / name
    folder.mkdir()
    input_path = folder / "carnivora.nwk"
    input_path.write_text((source.rstrip() + "\n") * multiplier)
    print(f"{name}: {125 * multiplier} genes, 53 taxa, {repeats} paired runs", flush=True)
    for repeat in range(repeats):
        prefixes = {}
        # Alternate order to reduce systematic warm-cache/order bias.
        for implementation in (implementations if repeat % 2 == 0 else list(reversed(implementations))):
            prefix = folder / f"{implementation.replace('+', 'p').replace(' ', '-')}-{repeat + 1}"
            prefixes[implementation] = prefix
            command = (
                ["Rscript", "--vanilla", str(repo / "tests/benchmark-worker.R"),
                 str(reference), str(input_path), str(prefix)]
                if implementation == "R"
                else [str(baseline if implementation == "C++ dense" else binary),
                      "run", "--trees", str(input_path), "--out", str(prefix), "--threads", "1"]
            )
            with Path(str(prefix) + ".stdout").open("w") as stdout, Path(str(prefix) + ".stderr").open("w") as stderr:
                subprocess.run(
                    ["/usr/bin/time", "-f", "%e\t%U\t%S\t%M", "-o", str(prefix) + ".time.tsv", *command],
                    env=environment, stdout=stdout, stderr=stderr, check=True, timeout=180,
                )
            elapsed, user, system, rss = Path(str(prefix) + ".time.tsv").read_text().split()
            analysis = (
                float(Path(str(prefix) + ".analysis-seconds.txt").read_text())
                if implementation == "R" else None
            )
            rows.append(dict(dataset=name, implementation=implementation, repeat=repeat + 1,
                             elapsed_seconds=float(elapsed), user_seconds=float(user),
                             system_seconds=float(system), peak_rss_kib=int(rss),
                             r_analysis_seconds=analysis))
            print(f"  {implementation}: elapsed={elapsed}s peak RSS={int(rss)/1024:.1f} MiB", flush=True)
        r_prefix = prefixes["R"]
        r_out = read_table(Path(str(r_prefix) + ".outliers.tsv"))
        r_scores = read_table(Path(str(r_prefix) + ".scores.tsv"))
        for implementation in implementations[1:]:
            cpp_prefix = prefixes[implementation]
            cpp_out = read_table(Path(str(cpp_prefix) + ".outliers.tsv"))
            assert r_out == cpp_out, f"{name}/{implementation}: outlier identities or order differ"
            cpp_scores = read_table(Path(str(cpp_prefix) + ".scores.tsv"))
            assert len(r_scores) == len(cpp_scores), f"{name}: state counts differ"
            assert [x["state"] for x in r_scores] == [x["state"] for x in cpp_scores]
            error = max(abs(float(r["quality"]) - float(c["quality"])) for r, c in zip(r_scores, cpp_scores))
            assert error < 1e-10, f"{name}/{implementation}: score mismatch {error}"
            print(f"  VERIFIED {implementation}: {len(r_out)} outliers, {len(r_scores)} states; max score error {error:.3g}", flush=True)
    item = dict(dataset=name, genes=125 * multiplier, taxa=53, repeats=repeats,
                outliers=len(r_out), states=len(r_scores))
    for implementation in implementations:
        samples = [r for r in rows if r["dataset"] == name and r["implementation"] == implementation]
        item[implementation] = {
            "median_elapsed_seconds": statistics.median(r["elapsed_seconds"] for r in samples),
            "elapsed_range_seconds": [min(r["elapsed_seconds"] for r in samples), max(r["elapsed_seconds"] for r in samples)],
            "median_peak_rss_mib": statistics.median(r["peak_rss_kib"] for r in samples) / 1024,
        }
        if implementation == "R":
            item[implementation]["median_analysis_seconds"] = statistics.median(r["r_analysis_seconds"] for r in samples)
    report.append(item)
    # Preserve completed datasets even if a subsequent scaling test fails.
    (destination / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    with (destination / "runs.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
print(json.dumps(report, indent=2))
