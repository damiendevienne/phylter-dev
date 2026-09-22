# Matrix-free RV benchmark — 2026-09-21

The default native CLI now obtains only the leading RV eigenpair, without
constructing RV. The older dense C++ executable was preserved before changing
the engine, and all three implementations were rerun in the same benchmark.

## Elapsed time and peak RAM

Medians from five sequential runs for Carnivora and three runs for each larger
case, one thread. Elapsed time includes startup, input reading, analysis and
compact output; peak RAM is resident process memory measured by GNU time.

| Genes | Original R | Previous dense C++ | Matrix-free C++ | R peak RAM | Dense C++ peak RAM | Matrix-free peak RAM |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 125 | 2.50 s | 0.28 s | **0.18 s** | 286.0 MiB | 16.8 MiB | **16.7 MiB** |
| 500 | 7.56 s | 3.77 s | **0.69 s** | 426.0 MiB | 49.4 MiB | **45.3 MiB** |
| 1,000 | 19.35 s | 25.18 s | **1.69 s** | 561.9 MiB | 108.6 MiB | **83.4 MiB** |

At 1,000 genes, the new solver makes the CLI approximately **14.9× faster than
the previous native version**, and **11.4× faster than a fresh original-R
process**. Peak RAM falls by approximately **23% relative to the previous
native version**.

R analysis time excluding startup/input loading was 1.150 s, 6.033 s and
17.318 s respectively. Much of the R-versus-C++ RAM difference reflects R's
interpreter, packages and larger retained result objects. The new-versus-old
C++ comparison isolates the RV solver change much more closely.

All inputs contain **53 taxa**. The 500/1,000-gene datasets repeat Carnivora
four/eight times: these are synthetic gene-count scaling tests, not independent
biological datasets or an evaluation of increasing taxon counts.

Observed elapsed ranges (R / previous C++ / matrix-free C++):

* 125 genes: 2.44–2.63 / 0.27–0.29 / 0.17–0.20 seconds.
* 500 genes: 7.47–7.59 / 3.68–3.87 / 0.69–0.71 seconds.
* 1,000 genes: 18.71–20.01 / 24.48–26.83 / 1.68–1.76 seconds.

No full diagnostic matrices were requested during the benchmark.
`--diagnostics` deliberately constructs RV for output and would conceal part
of the performance gain.

## Correctness

Every benchmark run matched the R reference and saved C++ baseline for exact
outlier identities/removal order, accepted-state counts and quality trajectories.
There were 94, 376 and 752 outliers, with 11 states in each dataset. Matrix-free
quality differences from R were at most 1.22e-15 in this benchmark.

Separately, all **18 R-reference configurations** passed, comparing WR,
weights, RV diagnostics, compromise and coordinate Gram matrices at every
accepted state. Maximum absolute array difference from R: **2.423e-12**.
Comparison with preserved dense C++ diagnostic outputs also passed all 18
cases; maximum absolute array difference: **2.480e-12**.
The comparison tolerance remains 1e-8 times max(1, reference magnitude).

The native test suite additionally checks explicit-versus-implicit products
(including asymmetric input), rank deficiency, negative RV entries, restarts,
close/tied leading eigenvalues, a leading vector orthogonal to the constant
start, explicit non-convergence and a 5,000-gene case without RV storage.
Serial/four-thread output equivalence and memory sanitizers are also checked.

## Reproduction

From the parent repository, with the unchanged original R installation in
`.audit/reference-lib` and the old C++ binary preserved before rebuilding:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tests/benchmark-resources.py \
  build/new-matrix-free-benchmark --multipliers 1 4 8 \
  --baseline build/baseline-dense/phylter
```

The recorded output directory is `build/resource-benchmark-matrix-free/`,
containing raw per-run logs, measurements, results, `runs.tsv`, `summary.json`
and SHA256 provenance in `environment.json`.
The saved dense baseline SHA256 is
`f28d90dbec4920be0ef57e2250546207f19ea9a88f82f91c623754e0e71d1362`.
The benchmark used the local Intel Core i7-1165G7 workstation, reference
BLAS/LAPACK and one thread; results are not hardware-independent guarantees.

See [RV.md](RV.md) for the operator, solver, convergence criterion and remaining
O(K N^2) memory requirement for distance/centered matrices.
