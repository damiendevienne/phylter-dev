# Phylter — standalone command-line application

Phylter finds gene/species pairs whose evolutionary signal is unusually
different from the rest of a phylogenomic dataset. It is a native C++
command-line successor to the original PhylteR R package and does not require
R, Python, or CRAN packages to run.

> **Project status:** Phylter 0.1.0 is available as a standalone download for
> Linux x86-64, macOS Apple Silicon, and Windows x86-64. This is the first
> development release; the original PhylteR R package remains the numerical
> reference. A Bioconda package is not yet available.

## Install a precompiled release

Download the archive for your system from the
[Phylter releases page](https://github.com/damiendevienne/phylter-dev/releases/latest):

| System | Archive |
| --- | --- |
| Linux, x86-64 | [`phylter-0.1.0-linux-x86_64.tar.gz`](https://github.com/damiendevienne/phylter-dev/releases/download/v0.1.0/phylter-0.1.0-linux-x86_64.tar.gz) |
| macOS, Apple Silicon (arm64) | [`phylter-0.1.0-macos-arm64.tar.gz`](https://github.com/damiendevienne/phylter-dev/releases/download/v0.1.0/phylter-0.1.0-macos-arm64.tar.gz) |
| Windows, x86-64 | [`phylter-0.1.0-windows-x86_64.zip`](https://github.com/damiendevienne/phylter-dev/releases/download/v0.1.0/phylter-0.1.0-windows-x86_64.zip) |

The archives contain Phylter and its required runtime libraries. They do not
require R, Python, Conda, BLAS, LAPACK, or a C++ compiler.

On Linux, download and extract the archive, then run:

```sh
tar -xzf phylter-0.1.0-linux-x86_64.tar.gz
./phylter-0.1.0-linux-x86_64/bin/phylter --version
```

On an Apple Silicon Mac:

```sh
tar -xzf phylter-0.1.0-macos-arm64.tar.gz
./phylter-0.1.0-macos-arm64/bin/phylter --version
```

The macOS binary is not currently signed or notarized. If macOS blocks it,
right-click the executable in Finder, select **Open**, and confirm that you
want to run it. Intel Macs are not yet supported by a precompiled archive;
follow the source-build instructions below instead.

On Windows, extract the ZIP archive and run this from PowerShell:

```powershell
.\phylter-0.1.0-windows-x86_64\bin\phylter.exe --version
```

Keep the archive's directory structure intact: the executable and bundled
libraries are designed to remain together. To invoke `phylter` from anywhere,
add its `bin` directory to your `PATH`.

Each archive has a matching `.sha256` file on the release page. After
downloading both files, Linux users can verify the download with
`sha256sum -c FILE.sha256`; on macOS use `shasum -a 256 -c FILE.sha256`.

## Build from source

Building is useful for development or for a platform without a precompiled
archive. You need:

- CMake 3.21 or newer
- a compiler with C++20 support (GCC, Clang, or recent MSVC)
- BLAS and LAPACK development libraries
- Git, if you clone the repository

On Ubuntu or Debian:

```sh
sudo apt update
sudo apt install build-essential cmake libblas-dev liblapack-dev
```

On macOS with Homebrew:

```sh
brew install cmake openblas lapack
```

Windows source builds are supported through CMake, but dependency setup varies.
The automated Windows build uses Conda; see
[the build workflow](.github/workflows/native.yml) for the exact dependencies.

Clone the repository and build an optimized executable:

```sh
git clone https://github.com/damiendevienne/phylter-dev.git
cd phylter-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable is now `build/phylter` on Linux and macOS.
Multi-configuration generators, including the usual Visual Studio setup, may
place it in `build/Release/` instead.

You can run it directly from the build directory, or install it under a local
prefix:

```sh
cmake --install build --prefix "$PWD/install" --component Runtime
export PATH="$PWD/install/bin:$PATH"
phylter --version
```

The `export` command affects only the current shell. Add the resulting `bin`
directory to your usual `PATH` configuration if you want a permanent install.
The installed executable uses the BLAS/LAPACK libraries available on the
system.

## Test a source build

Run the complete test suite from the repository root:

```sh
ctest --test-dir build --output-on-failure
```

There are five tests covering the core algorithms, frozen R reference results,
input validation, command-line output, and parallel reproducibility. R is not
required. A successful run ends with:

```text
100% tests passed, 0 tests failed out of 5
```

## Check your installation and input data

Every release archive contains the Carnivora example dataset. From inside the
extracted archive, check both the executable and this input with:

```sh
bin/phylter --check --trees share/phylter/examples/carnivora.nwk
```

On Windows, use `bin\phylter.exe` and Windows-style paths instead. A successful
check ends with `Input is valid and ready for analysis.`

From a source build, the equivalent command is:

```sh
./build/phylter --check --trees tests/fixtures/carnivora.nwk
```

This verifies the Newick structure, labels, and supplied branch lengths. It
also reports whether patristic distances and support-based branch collapsing
are available, and shows the gene identifiers derived from the filenames.

`--check` accepts only `--trees`; it does not predict whether a particular set
of optimization parameters will produce outliers.

## Run the included example

The release archive includes a real Carnivora dataset. From its extracted
directory, run:

```sh
mkdir -p results
bin/phylter \
  --trees share/phylter/examples/carnivora.nwk \
  --out results/carnivora
```

In PowerShell, use `New-Item -ItemType Directory -Force results` and replace
`bin/phylter` with `bin\phylter.exe`. For a source build, the corresponding
input is `tests/fixtures/carnivora.nwk`.

During the run, Phylter reports the dataset dimensions and taxon coverage,
branch-length and node-support availability, active settings, and every cell or
complete-gene proposal with its score change and acceptance decision. The final
summary reports quality gain, filtered-data proportion, complete outliers,
elapsed time, peak memory, and the exact paths of the four result files. Before
optimization, a rough upper memory estimate is compared with currently
available RAM and a warning is shown for potentially unsafe runs:

- `results/carnivora.outliers.tsv` — outlier gene/species pairs, in removal order
- `results/carnivora.discarded.tsv` — pairs excluded by the normalization cutoff
- `results/carnivora.scores.tsv` — quality score for every accepted state
- `results/carnivora.summary.txt` — result counts and the parameters used

The default analysis finds 94 outlier gene/species pairs in 10 accepted
iterations, with a final quality of approximately `0.944259987303`.

## Analyze your own trees

Pass either one Newick file or a directory containing Newick files:

```sh
phylter --trees gene_trees/ --out results/analysis --threads 4
```

The output prefix can include a directory; Phylter creates it when necessary.
It will not overwrite existing result files unless you add `--force`.

Supported filename extensions are `.nwk`, `.newick`, `.tre`, `.tree`,
and `.treefile` (case-insensitive). A file may contain one or several trees.
When a directory is supplied, files are read in lexicographic filename order.

By default, Phylter uses patristic distances, which require finite,
non-negative branch lengths. For trees without branch lengths, use nodal
distances:

```sh
phylter \
  --trees gene_trees/ \
  --distance nodal \
  --out results/nodal-analysis
```

Quoted labels, bracket comments, multifurcations, and datasets in which genes
contain different sets of taxa are supported. Each input gene must contain at
least three taxa, and at least two genes must remain after filtering.

Gene identifiers come from filenames. A file named `gene42.nwk` produces
`gene42`. If it contains several trees, their identifiers are `gene42:1`,
`gene42:2`, and so on. Taxon labels, including underscores, are preserved.
Identifiers cannot contain tabs or newlines.

To inspect tree structure and see which analyses the dataset supports:

```sh
phylter --check --trees gene_trees/
```

## Use precomputed distance matrices

Instead of trees, Phylter can read a directory containing one square TSV
distance matrix per gene:

```sh
phylter --matrices distance_matrices/ --out results/analysis
```

Each `.tsv` file must have `taxon` in the first header cell, followed by the
taxon names. Every data row starts with its corresponding taxon name. Row and
column order must match. Values must be finite, symmetric, and non-negative,
with zeros on the diagonal.

```text
taxon	A	B	C
A	0	0.5	0.8
B	0.5	0	0.4
C	0.8	0.4	0
```

Matrices may contain different sets of taxa; Phylter imputes missing values
before the DISTATIS analysis.

## Common options

Run `phylter --help` for the complete list.

| Option | Purpose |
| --- | --- |
| `--threads N` | Use `N` workers for RV computations (default: 1) |
| `--distance patristic\|nodal` | Select the distance calculated from trees |
| `--norm median\|mean\|none` | Select distance normalization |
| `--support-cutoff N` | Collapse internal branches below this node-support threshold |
| `--initial-only` | Calculate only the initial state |
| `--force` | Overwrite existing result files |
| `--diagnostics DIR` | Save each state's matrices; uses substantial disk space |

Result files report outliers; Phylter does not modify or prune the original
trees or alignments.

## Troubleshooting

**CMake cannot find BLAS or LAPACK.** Install their development packages, not
only the runtime libraries. On Debian/Ubuntu these are `libblas-dev` and
`liblapack-dev`.

**An output file already exists.** Choose another `--out` prefix, remove the
old results, or use `--force` if overwriting is intentional.

**`--threads` reports that OpenMP is unavailable.** Use `--threads 1`, or
configure again with an OpenMP-capable toolchain. A threaded BLAS may create its
own workers; use `OPENBLAS_NUM_THREADS` or `MKL_NUM_THREADS` to control the
total thread count.

**You moved the source directory after configuring CMake.** Delete the build
directory and rerun the configure and build commands. CMake build directories
are not portable between source locations.

## Packaging and validation details

Release archives are built and smoke-tested automatically on Linux, macOS, and
Windows before publication. A development Conda recipe is also included, but it
has not been published to Bioconda. See
[packaging/README.md](packaging/README.md) for maintainer details.

The regression fixtures were exported from R reference revision
`4d74241169be3882da9d5f08da47bd40604764eb`. The port preserves the reference
algorithm's normalization, DISTATIS weighting, eigenvalue selection, outlier
thresholds, clustering, and proposal acceptance behavior. Comparisons are
tolerance-based because LAPACK implementations can differ slightly.

Developers can find the full numerical methodology in
[VALIDATION.md](VALIDATION.md). Performance and implementation notes are in:

- [RV solver design](RV.md)
- [matrix-free solver benchmark](BENCHMARK-MATRIX-FREE.md)
- [real Carnivora performance and memory benchmark](BENCHMARK-REAL-CARNIVORA.md)
- [gene and taxon scaling benchmark](BENCHMARK-TIP-SCALING.md)

The application still stores dense distance and centered matrices. Packed
symmetric storage and an approximate fast mode are not implemented. Diagnostic
output additionally materializes the full RV matrix, so it uses more memory
than a normal run.
