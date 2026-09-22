# C++ versus the original R implementation — 2026-09-21

Historical benchmark of the first, dense-RV C++ port. The default solver has
since changed to [matrix-free RV](RV.md); these C++ timings describe the saved
previous executable, not the new default.
See the [updated comparison](BENCHMARK-MATRIX-FREE.md) for the new measurements.

The native CLI uses substantially less RAM in all tested cases. It is faster
on Carnivora and its 500-gene replication, but **slower at 1,000 genes**.
The first port's speed advantage does not extend uniformly to larger inputs.

## Measurements

Median elapsed times and median peak resident RAM, measured in separate fresh
processes using GNU time. Runs are sequential, alternating R/C++ order, with
one thread and the same numerical settings. R uses `parallel=FALSE`; the CLI
uses `--threads 1`. BLAS thread environment variables are also set to one.

| Input | Genes | Original R, elapsed | C++, elapsed | Original R, peak RAM | C++, peak RAM |
| --- | ---: | ---: | ---: | ---: | ---: |
| Carnivora | 125 | 2.38 s | 0.28 s | 285.9 MiB | 16.8 MiB |
| Carnivora repeated 4 times | 500 | 7.40 s | 3.79 s | 424.2 MiB | 49.5 MiB |
| Carnivora repeated 8 times | 1,000 | 17.34 s | 23.31 s | 562.2 MiB | 108.5 MiB |

All inputs contain **53 taxa**. The larger inputs duplicate the same trees;
they are synthetic tests of gene-count scaling, not independent biological
datasets or tests of increasing taxon count.

Five paired runs were used for 125 genes and three paired runs for each larger
input. Observed elapsed ranges, R / C++:

* 125 genes: 2.36–2.39 s / 0.27–0.29 s.
* 500 genes: 7.35–7.46 s / 3.73–3.80 s.
* 1,000 genes: 16.79–17.83 s / 22.49–23.35 s.

The end-to-end speed ratios are approximately **8.5× faster**, **2.0× faster**,
and **1.34× slower** for C++, respectively. Its peak RAM is approximately
17.0×, 8.6× and 5.2× lower.

## Startup versus analysis time

The main table includes process startup, package loading, input reading,
analysis and compact result writing. The original R function also retains its
full initial/final result objects; C++ uses its compact default state.

The R analysis alone, measured after package loading and tree reading, took
median times of **1.124 s**, **5.899 s**, and **15.477 s**. The C++ times above
still include its input/output and startup. Thus the original Carnivora
analysis speedup is approximately **4×** when R startup is excluded, consistent
with the earlier benchmark.

Peak R RAM includes its interpreter and loaded packages. The memory reduction
does not arise solely from more efficient matrix storage: both implementations
still use dense matrices, and fixed R overhead matters greatly on small inputs.

## Correctness and interpretation

Every paired run checked exact outlier identities **and removal order**, the
number/order of accepted states, and the entire quality trajectory:

| Genes | Outlier pairs | Accepted states | Largest absolute score difference |
| ---: | ---: | ---: | ---: |
| 125 | 94 | 11 | 1.11e-15 |
| 500 | 376 | 11 | 2.11e-15 |
| 1,000 | 752 | 11 | 2.66e-15 |

The original reference is PhylteR 0.9.12 at revision
`4d74241169be3882da9d5f08da47bd40604764eb`, installed independently in
`.audit/reference-lib`. The scientific engine was not modified during this
benchmark.

A likely reason for the crossover is visible in the source: the C++ port
currently calls a full RV eigendecomposition (`symmetric_eigen(rv)`, LAPACK
DSYEV), while the original R code requests only the leading eigenpair using
`eigs_sym(RVmat, 1)`. A full decomposition has cubic cost in gene count.
This is a code-based explanation, not a timing breakdown from a profiler.
Computing only the leading eigenpair, then introducing a matrix-free RV
operator, is a priority before claiming improved speed on large datasets.

## Reproduce

On the tested Linux workstation (Intel Core i7-1165G7), from the parent
repository, with the reference R installation and release C++ binary available:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tests/benchmark-resources.py \
  build/new-resource-benchmark

PYTHONDONTWRITEBYTECODE=1 python3 tests/benchmark-resources.py \
  build/new-resource-benchmark-1000 --multipliers 8 --repeats 3
```

Use fresh output directories. Scripts write only inside the repository.
Python and R are used for this developer benchmark, not by the native CLI.
Each directory contains per-run logs, GNU time measurements, score/outlier
files, `runs.tsv` and `summary.json`. The recorded runs are in
`build/resource-benchmark-20260921/` and
`build/resource-benchmark-1000-20260921/`; their `environment.json` files
record the native binary and input SHA256 hashes.
