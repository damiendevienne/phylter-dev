# Distribution status

The Linux executable and bundled archive can be built and tested locally.
The Conda recipe is a development recipe, not a published Bioconda package.
macOS and Windows builds still need execution on those operating systems.
Do not describe unexecuted builds as validated releases.

## Linux archive

From the standalone project's root:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DPHYLTER2_BUNDLE_RUNTIME=ON
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
cpack --config build-release/CPackConfig.cmake -B build-release
```

Keep the extracted bin/ and lib/ directories together. Run bin/phylter directly.
The test data are included in share/phylter/examples/carnivora.nwk:

```sh
bin/phylter run --trees share/phylter/examples/carnivora.nwk --out results/carnivora
```

The archive includes numerical/compiler runtimes, but uses the host's glibc;
build public Linux releases on the oldest supported Linux baseline.

Before redistribution, collect license notices and corresponding source
information for the bundled runtime versions. The local archive is a development
artifact, not yet a public release.

## Conda / future Bioconda

```sh
conda build packaging/conda -c conda-forge
```

Only the build environment needs CMake, Ninja and a compiler. The resulting
runtime package needs numerical libraries, not R or Python. OpenMP is disabled
in this first recipe to avoid an untested compiler/runtime combination.
When using this inside the original R repository, set Conda's cache, environment
and build directories inside the repository if following the local-only policy.

Once the dedicated repository exists, use an immutable source archive with its
SHA256 and add the project's homepage, documentation and maintainers to the
recipe before submitting it. Bioconda publication and cross-platform archives
are separate release tasks.
