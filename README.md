# PhylteR 2 — native C++ CLI

The complete filtering loop now runs as a native executable, with no R, Python
or CRAN packages at runtime. This is a development version of the standalone
port. The numerical reference remains the original PhylteR R package.

## Try the included dataset

From this standalone project's directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/phylter run --trees tests/fixtures/carnivora.nwk --out build/carnivora
```

From the **parent R repository**, use:

```sh
./build/phylter run \
  --trees tests/fixtures/carnivora.nwk \
  --out build/carnivora
```

Expected default result: **94 outlier gene/species pairs**, 10 accepted
iterations (11 states), final quality **0.944259987303**.

Install locally and put the executable on your PATH:

```sh
cmake --install build --prefix "$PWD/install" --component Runtime
export PATH="$PWD/install/bin:$PATH"
phylter run --trees gene_trees/ --out results/analysis --threads 4
```

Building requires CMake >=3.21, a C++20 compiler and LP64 BLAS/LAPACK development
libraries. A prebuilt Linux bundle includes numerical libraries: extract it,
keep bin/ and lib/ together, and run bin/phylter.
See [packaging/README.md](packaging/README.md). No package is published on
Bioconda yet; macOS/Windows executables have not yet been validated.

## Input and output

Use `phylter --help` for all options. A Newick file may contain several trees;
directories are read in lexicographic filename order. Supported extensions
(case insensitive): .nwk, .newick, .tre, .tree, .treefile.
Quoted labels, bracket comments and multifurcations are accepted.
Patristic distances require nonnegative finite branch lengths. Nodal distances
also accept trees without branch lengths.

Gene identifiers are the filename stem for one tree per file, or
`stem:1`, `stem:2`, etc. for multiple trees. Taxon identifiers are preserved,
including underscores. Identifiers cannot contain tabs or newlines.
This naming convention differs from the older R CLI's numeric IDs for a
multi-tree file; biological results are compared using an explicit gene mapping.

You may provide matrices using `--matrices DIRECTORY`: one square matrix per
.tsv file, first header cell `taxon`, then column taxon names, and each data row
starting with the corresponding row name. Row/column orders must match.
Distances must be finite, symmetric, nonnegative, with zero diagonals.
Genes may contain different taxa: imputation makes their dimensions and taxon
order identical **before** DISTATIS. At least two retained genes and three taxa
per input gene are required.

Compact results:

* `PREFIX.outliers.tsv`: gene/species pairs in removal order.
* `PREFIX.discarded.tsv`: pairs discarded by the normalization cutoff.
* `PREFIX.scores.tsv`: quality at every accepted state.
* `PREFIX.summary.txt`: counts, complete outliers and analysis parameters.

These are outlier reports; the CLI does not prune original trees or alignments.
Files are protected against overwriting unless `--force` is supplied.
`--diagnostics DIRECTORY` additionally writes WR, RV, compromise, weights and
coordinate Gram matrices at every accepted state. This can consume substantial
disk space. The diagnostics directory must be empty.

## Numerical compatibility

The port preserves distance normalization, double centering of **unsquared**
distances, RV weighting, largest-magnitude eigenvalue selection, broken-stick
axis selection, WR normalization, adjusted Tukey/medcouple thresholds,
complete-link clustering and its tie order, island maxima, the whole-gene pass,
and acceptance/rejection of proposed removals.

Two legacy behaviors are preserved deliberately:

* Never-co-occurring taxa use R's existing neighbor fallback. The original
  `mean(a,b)` expression returns scalar `a` because `b` is the trim argument;
  the resulting imputed matrices can be asymmetric. Disconnected pairs with no
  suitable shared neighbor fail explicitly. A revised symmetric estimator needs
  a separately validated method/version.
* `--support-cutoff` reproduces the existing R rule: low-support **parent**
  nodes modify outgoing edges, then short internal edges collapse. This is not
  a new interpretation of bootstrap filtering.

Numerical equality is tolerance-based; LAPACK and RSpectra need not return
identical floating-point bits or eigenvector signs. Near a threshold, or across
a tied eigenvalue cutoff, platform-dependent numerical differences remain
possible. The default regression compares exact outlier IDs/order and scores,
WR, RV, weights, compromise and coordinate geometry.

## Validation

`ctest` uses frozen fixtures exported from the unchanged R reference revision
`4d74241169be3882da9d5f08da47bd40604764eb`. **R is not needed for these tests.**
It also checks medcouple samples, clustering ties, invalid inputs, matrix size
overflow, Newick parsing, CLI output and serial/parallel reproducibility.

Developers with the original R package installed can rerun the larger comparison:

```sh
# From the parent R repository; all temporary files remain in that repository.
TMPDIR="$PWD/.audit/tmp" \
R_LIBS="$PWD/.audit/reference-lib:$PWD/.audit/library" \
TZ=Europe/Athens Rscript tests/compare-reference.R \
  build/phylter build/new-reference-validation
```

Use a fresh output directory. This compares **every accepted state** in 18
configurations, including missing taxa, never-co-occurring taxa, negative
spectra, support filtering, rejected proposals and no-outlier cases.
See [VALIDATION.md](VALIDATION.md), the
[matrix-free solver benchmark](BENCHMARK-MATRIX-FREE.md), and the
[real Carnivora speed/RAM comparison and laptop estimate](BENCHMARK-REAL-CARNIVORA.md).
The [genes × tips benchmark](BENCHMARK-TIP-SCALING.md) expands the input trees
to 106/212 taxa and measures both dimensions, with reproducible synthetic inputs.
The default RV solver is now [matrix-free](RV.md); `--rv-method dense` retains
the previous solver for comparison. Diagnostic output still requests a full RV
matrix and therefore has a higher memory/time cost than a normal run.

## Performance and next stages

Already implemented:

* contiguous matrix storage and BLAS matrix products;
* no projections/eigendecomposition of the compromise for rejected candidates;
* compact default state retention (no history of large matrices in memory);
* one-gene-at-a-time expansion during imputation;
* deterministic OpenMP parallelism over independent entries of RV products
  (or RV pairs with the dense backend);
* the leading RV eigenpair from an implicit operator, with no full RV allocation
  or eigendecomposition in the default CLI path.

`--threads` controls the RV workers (default 1). BLAS calls occur outside the
OpenMP region, avoiding nested parallelism. A threaded BLAS has its own thread
setting; use `OPENBLAS_NUM_THREADS` or `MKL_NUM_THREADS` to match the intended
budget. The tested Linux build uses the sequential reference BLAS.

The port still stores dense distance/centered matrices. Packed symmetric
storage and an explicit approximate `fast` mode are **not implemented yet**. They are the next
separately benchmarked stages, using this port and R as references. The
[real-data benchmark](BENCHMARK-REAL-CARNIVORA.md) measures up to 3,000 genes
and estimates full-data memory at 53 taxa; larger taxon counts need separate tests.

ERABLE-based alternatives, incremental updates and GPU computation remain
outside this stage. Plots are intentionally absent.
