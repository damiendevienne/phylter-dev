# Real Carnivora benchmark — 2026-09-21

On the same random sample of 1,000 real genes, the native matrix-free CLI is
**9.33× faster** than original R and uses **86.4% less peak resident memory**.
All reported outlier identities, their order, discarded entries and accepted
quality trajectories agree. No scientific-engine change was needed for this test.

## Data and reference version

Source: [Carnivora trees archive](https://github.com/damiendevienne/phylter-data/blob/main/Carnivora/data/trees.tar.xz).
The downloaded archive contains **14,463 genes**, **53 taxa overall**, and
5–53 taxa per gene (median 47). Unlike the earlier repeated-gene benchmarks,
these samples contain distinct genes selected without replacement.

The archive, samples and per-gene manifests are stored under
`validation-data/carnivora-full/`, ignored by Git. The archive is retained without
extracting thousands of files. Sampling preserves the original Newick text and
branch-length precision. Seed **20260921** selects nested samples from a shuffled
list of sorted archive paths; selected trees are then ordered by archive path.
Manifests map identifiers such as `sample-1000:242` to the original gene filename.
The census only reads tree labels; it does not construct full-data distance matrices.

Archive SHA256:
`be9c03f2c2e2a4a303302d5a834cc6071ac2cc8fd2371c4e8c7af7da3f95a694`.

1,000-gene input SHA256:
`0694f26427c51422af7a1aaebd998c68ee5b935b53d8d078fb67524e4854d146`.

**R is the unchanged original version, not the optimized R branch.** This
repository calls its primary branch `master`, rather than `main`. Local `master`
and `origin/master` both identify commit
`4d74241169be3882da9d5f08da47bd40604764eb`, package version 0.9.12.
The archived reference's 19 R source files were checked byte-for-byte against
that commit. Each worker explicitly loads `.audit/reference-lib/phylter`, and
records the loaded package path/version. The optimized R package in
`.audit/candidate-lib` is not used in these measurements.

Native executable SHA256:
`8d6b994899fc937e85ca530433f68f17e80f157994f190c16799c6e6dcf347d5`.

## Comparable work and measurements

Both implementations read exactly the same multi-tree Newick file and perform
the complete filtering analysis. Parameters: patristic distances, median
normalization, cutoff 0.001, k = k2 = 3, row WR normalization, islands enabled,
minimum quality gain 1e-5, support cutoff 0. Missing taxa are imputed before
DISTATIS, producing the same **53 × 53** dimensions for every gene.

Runs use one thread (`parallel=FALSE` in R, `--threads 1` in C++, BLAS thread
environment variables fixed to 1), the same machine and reference BLAS.
Hardware: Intel Core i7-1165G7, nominal 16 GB RAM (14.86 GiB visible to Linux).
The C++ executable is a Release build using system BLAS/LAPACK.

GNU time measures elapsed time and maximum resident set size of a fresh process.
Elapsed time includes startup, reading trees, analysis and compact result output.
Runs are sequential; R/C++ order alternates across repetitions. Neither version
plots. C++ uses default compact output and matrix-free RV, without diagnostics.
R still retains its normal rich result object internally, including RV and
projection matrices: this is a comparison of complete workflows producing the
same filtering results, not an isolated measurement of programming languages.

| Genes | Repetitions per implementation | Original R | C++ | R peak RAM | C++ peak RAM |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 1, pilot | 2.90 s | 0.22 s | 282.4 MiB | 14.7 MiB |
| 250 | 1, pilot | 5.84 s | 0.73 s | 332.5 MiB | 25.7 MiB |
| **1,000** | **3** | **43.85 s** | **4.70 s** | **597.7 MiB** | **81.2 MiB** |

Three-run figures are medians. At 1,000 genes, elapsed ranges were
43.15–43.89 seconds for R and 4.68–4.71 seconds for C++.
Maximum observed peak RSS was 597.70 MiB for R and 81.29 MiB for C++.
R analysis alone, excluding package startup and Newick reading, took a median
**41.376 seconds**, so startup accounts for only a small part of the difference.

Every process had a **3 GiB address-space limit** and a 600-second timeout.
The runner required 5 GiB available before starting and monitored available
system memory, stopping the process group if it dropped below 2 GiB.
Core dumps were disabled. All runs completed normally; no limit was reached.
All temporary files, logs and outputs stayed within this repository.

## Correctness

Every paired run passed these checks:

* Exact outlier gene/species pairs **and removal order**.
* Exact discarded gene/species entries (none in these samples).
* Identical accepted-state count and identifiers.
* Every accepted quality score within 1e-10; the largest observed difference
  was **1.11e-15**.

| Genes | Outlier pairs | Accepted states, including initial |
| ---: | ---: | ---: |
| 100 | 92 | 9 |
| 250 | 171 | 15 |
| 1,000 | 538 | 23 |

For 1,000 genes, initial quality was 0.792569818554466 and final quality was
0.849392933229086. The 23 states represent the initial state plus 22 accepted
iterations. This benchmark compares compact scientific results and trajectories;
it does not export/compare every intermediate matrix. See
[VALIDATION.md](VALIDATION.md) for the separate matrix-level regression suite.

## Full dataset: laptop feasibility without running it

**No full 14,463-gene analysis was launched in R or C++.**

An additional C++-only scaling check used **3,000 genes**, still 53 taxa, under
the same memory limits. Three runs took 14.35–14.70 seconds (median **14.36 s**)
and a median **225.94 MiB** peak RSS (maximum 226.14 MiB). All three produced
1,344 outlier pairs and 24 accepted states. These runs are a memory/scaling
check, **not a correctness comparison against R at 3,000 genes**.

At fixed taxon count, the default C++ implementation stores O(K × N²) distance
and centered matrices, plus smaller working arrays. It does not store the
K × K RV matrix. The code retains approximately three full matrix stacks
during an accepted proposal: current matrices, proposed matrices and centered
proposed matrices. Projection arrays, labels, WR and allocator overhead add to this.

For K = 14,463 and N = 53:

* One full stack of double-precision distance matrices is **0.303 GiB**.
* Three such stacks are **0.908 GiB** before other allocations.
* Scaling measured peak RSS directly from 1,000 and 3,000 genes gives about
  **1.15 GiB** and **1.06 GiB**, respectively.

Taken together, these support an approximate **1.1–1.3 GiB process RAM estimate**.
Allow **2 GiB for phylter** as a practical margin and start with at least
**3 GiB available system RAM**. A 16 GB laptop such as this one should have ample
capacity; an 8 GB laptop should also be feasible if enough memory is free.
This is an extrapolation supported by sample measurements and code inspection,
not a measured full-run maximum or a formal memory bound. A hard process limit
can stop a run if the estimate is exceeded.

The 3,000-gene elapsed time scales to about 69 seconds at 14,463 genes **if the
iteration count and solver convergence stay similar**. Actual time may differ;
the full analysis could require more iterations. Memory is the more defensible
feasibility estimate here.

These estimates apply specifically to **53 taxa**, the **default matrix-free
solver**, **compact output** and **one thread**. Do not extrapolate them to much
larger taxon sets: distance matrix storage grows quadratically with taxon count.
`--rv-method dense` or `--diagnostics` reintroduces full RV allocation: one such
matrix alone is **1.56 GiB** at 14,463 genes, before copies/workspaces. Original R
also allocates RV and several large temporaries; no full-data R run was attempted.

## Reproduce

From the parent R repository, using the existing original R reference installation
and built native CLI, choose fresh output directories:

```sh
python3 tests/prepare-carnivora.py
python3 tests/benchmark-carnivora.py \
  build/new-carnivora-comparison --sizes 1000 --repeats 3
python3 tests/benchmark-carnivora.py \
  build/new-carnivora-scaling --sizes 3000 --repeats 3 --cpp-only
```

Raw measurements, logs, results and provenance from this session are in:

* `build/carnivora-pilot/`
* `build/carnivora-benchmark/`
* `build/carnivora-scaling/`

Each contains `runs.tsv`, `summary.json` and `environment.json`. Data preparation
uses Python 3.14.4 here; Python/R are developer validation tools, not dependencies
of the native executable. To run the sampled analysis directly:

```sh
build/phylter run \
  --trees validation-data/carnivora-full/sample-1000.nwk \
  --out build/my-carnivora-result \
  --threads 1
```
