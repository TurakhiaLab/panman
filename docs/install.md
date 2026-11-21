# Installation on macOS

Install Apple's Command Line Tools if needed:

```bash
xcode-select --install
```

From a PanMAN checkout, build the current source (including local changes):

```bash
bash install/installationMacOS.sh --threads 2
./build/bin/panmanUtils --help
```

The installer supports Intel and Apple Silicon on macOS 11 or newer. It downloads
and verifies Miniforge, creates `build/source-env` with build dependencies,
and builds `panmanUtils` in the checkout’s `build/` directory. It uses classic TBB 2020.2
because the source still uses APIs removed in oneTBB. The wrapper at
`./build/bin/panmanUtils` runs the executable in its dependency environment.
No shell activation or sudo is required. Rerun the installer after source changes.

Use `--prefix /path/to/install` to change the installation directory (no spaces),
`--repo /path/to/panman` to select another checkout, and `--threads N` to control
build parallelism. The default checkout is resolved relative to the installer.

To install the published Bioconda package instead:

```bash
bash install/installationMacOS.sh --package
# Or choose a particular release:
bash install/installationMacOS.sh --package --version 0.1.4
```

The package uses a separate `build/env` directory by default. Both modes update the same command
wrapper, so the most recently completed installation selects the executable.
The published 0.1.4 package does not support MAF input (`-A`); use a source build
for that feature. Check `panmanUtils --help` for `--input-maf`.

# Installation on Ubuntu

From a PanMAN checkout on Ubuntu 22.04 or 24.04:

```bash
bash install/installationUbuntu.sh
./build/bin/panmanUtils --help
```

The installer supports x86_64 and aarch64. It installs missing bootstrap tools
(curl, CA certificates, bzip2, and build-essential) using apt-get, with sudo only
when needed. It verifies a pinned Miniforge installer and creates private build
dependencies in `build/source-env`. It builds the current checkout, including
local changes, in the checkout’s `build/` directory; it does not install libraries globally
or change shell initialization files. Run it as your normal user, rather than
running the entire script with sudo.

Use `--prefix DIR`, `--repo DIR`, or `--threads N` to customize the installation.
The prefix must not contain whitespace. Rerun the same command to rebuild after
source changes. `PREFIX/explicit-spec.txt` records the resolved dependencies.

Validate MAF import/export, including gap handling, with:

```bash
python3 tests/test_maf_export.py ./build/bin/panmanUtils
```

Ubuntu builds default to one compiler job and omit debug information to reduce
peak memory use. Increase `--threads` only when sufficient RAM is available.

For both source installers, the build directory is always `REPO/build`, where
`REPO` is the checkout containing the installer, or the checkout selected with
`--repo`. CMake creates this directory if needed and reuses it on subsequent
runs. The executable is `REPO/build/panmanUtils`. The `--prefix` option controls
Miniforge, dependencies, and the wrapper, and does not change the build location.
Use `PREFIX/bin/panmanUtils` to run with the dependency environment.

The default prefix is also `REPO/build`, resolved after `--repo` is processed.
By default, Miniforge is in `build/miniforge`, source dependencies are in
`build/source-env`, and the command wrapper is `build/bin/panmanUtils`.
An explicit `--prefix DIR` overrides these dependency and wrapper locations.
