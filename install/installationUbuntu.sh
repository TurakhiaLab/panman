#!/usr/bin/env bash
# Build the local PanMAN checkout with isolated, compatible dependencies on Ubuntu.
set -euo pipefail

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
usage() {
  cat <<'HELP'
Usage: bash installationUbuntu.sh [--prefix DIR] [--repo DIR] [--threads N]

Build the local PanMAN checkout on Ubuntu (x86_64 or aarch64).
Default prefix: REPO/build (must not contain spaces).
Source builds go in REPO/build; --prefix controls dependencies and the wrapper.
Only missing system bootstrap tools require root or sudo. Build dependencies
live in a private Miniforge environment; no shell activation is required.
Run PREFIX/bin/panmanUtils after installation. Rerun to rebuild source changes.
The default is one build job; increase --threads on machines with more RAM.
HELP
}
prefix=''
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
threads=1
while (($#)); do
  case "$1" in
    --prefix|--repo|--threads)
      (($# >= 2)) || fail "$1 requires a value"
      [[ -n $2 && $2 != --* ]] || fail "$1 requires a value"
      case "$1" in
        --prefix) prefix=$2 ;;
        --repo) repo=$2 ;; --threads) threads=$2 ;;
      esac
      shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown option: $1" ;;
  esac
done
[[ $threads =~ ^[1-9][0-9]*$ ]] || fail '--threads must be a positive integer.'
[[ $(uname -s) == Linux ]] || fail 'This installer requires Ubuntu Linux.'
[[ -f /etc/os-release ]] || fail 'Cannot identify the Linux distribution.'
. /etc/os-release
[[ ${ID:-} == ubuntu ]] || fail 'This installer targets Ubuntu.'
case "$(uname -m)" in
  x86_64) arch=x86_64; platform=linux-64 ;;
  aarch64|arm64) arch=aarch64; platform=linux-aarch64 ;;
  *) fail "Unsupported architecture: $(uname -m)" ;;
esac
[[ -f $repo/CMakeLists.txt && -f $repo/src/panmanUtils.cpp ]] || fail "Not a PanMAN source checkout: $repo"
repo=$(cd "$repo" && pwd -P)
prefix=${prefix:-"$repo/build"}
[[ $prefix != *[[:space:]]* ]] || fail 'Choose a prefix without whitespace (Miniforge requirement).'
# Install only bootstrap tools globally. Libraries and code generators are private.
missing=()
command -v curl >/dev/null || missing+=(curl)
[[ -f /etc/ssl/certs/ca-certificates.crt ]] || missing+=(ca-certificates)
command -v bzip2 >/dev/null || missing+=(bzip2)
command -v g++ >/dev/null && command -v make >/dev/null || missing+=(build-essential)
if ((${#missing[@]})); then
  elevate=()
  if ((EUID != 0)); then
    command -v sudo >/dev/null || fail "Install these system packages as root, then retry: ${missing[*]}"
    elevate=(sudo)
  fi
  "${elevate[@]}" apt-get update
  "${elevate[@]}" env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends "${missing[@]}"
fi
mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd -P)
[[ $prefix != *[[:space:]]* ]] || fail 'The resolved prefix contains whitespace.'
forge="$prefix/miniforge"
env_dir="$prefix/source-env"
# Ignore an active Conda environment or a platform override in the calling shell.
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SUBDIR PYTHONHOME PYTHONPATH
export CONDA_SUBDIR="$platform"
if [[ ! -x "$forge/bin/conda" ]]; then
  [[ ! -e $forge ]] || fail "Incomplete Miniforge installation at $forge; use a new prefix."
  download_dir=$(mktemp -d "${TMPDIR:-/tmp}/panman-install.XXXXXX")
  trap 'rm -rf "$download_dir"' EXIT
  # Checksum assets exist only for versioned installer filenames.
  miniforge_version=26.9.4-2
  installer="Miniforge3-${miniforge_version}-Linux-$arch.sh"
  url="https://github.com/conda-forge/miniforge/releases/download/$miniforge_version"
  curl -fL --retry 3 "$url/$installer" -o "$download_dir/$installer"
  curl -fL --retry 3 "$url/$installer.sha256" -o "$download_dir/checksum"
  expected=$(awk 'NR == 1 {print $1}' "$download_dir/checksum")
  [[ $expected =~ ^[a-fA-F0-9]{64}$ ]] || fail 'Invalid Miniforge checksum.'
  (cd "$download_dir" && printf '%s  %s\n' "$expected" "$installer" | sha256sum -c -)
  bash "$download_dir/$installer" -b -p "$forge"
fi
operation=create
if [[ -f $env_dir/conda-meta/history ]]; then
  operation=install
elif [[ -e $env_dir ]]; then
  fail "Not a Conda environment: $env_dir; use a new prefix."
fi
# PanMAN still uses task_scheduler_init and other pre-oneTBB APIs.
"$forge/bin/conda" "$operation" -y -p "$env_dir" \
  --override-channels --strict-channel-priority -c conda-forge \
  'cmake>=3.20,<4' make pkg-config 'libboost-devel>=1.84,<1.86' \
  capnproto=1.0.2 'protobuf>=3.20,<3.21' jsoncpp \
  tbb=2020.2 tbb-devel=2020.2 zlib xz
# Avoid the project's unconditional debug information in release builds.
# Link the environment's C++ runtime: packaged libraries may require a newer
# libstdc++ ABI than the system compiler provides (notably on Ubuntu 22.04).
build="$repo/build"
"$forge/bin/conda" run --no-capture-output -p "$env_dir" \
  cmake -S "$repo" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG -g0" \
  -DCMAKE_EXE_LINKER_FLAGS="-L$env_dir/lib" -DCMAKE_BUILD_RPATH="$env_dir/lib" \
  -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_PREFIX_PATH="$env_dir" \
  -DCAPNP_EXECUTABLE="$env_dir/bin/capnp" \
  -DCAPNPC_CXX_EXECUTABLE="$env_dir/bin/capnpc-c++"
"$forge/bin/conda" run --no-capture-output -p "$env_dir" \
  cmake --build "$build" --target panmanUtils --parallel "$threads"
executable="$build/panmanUtils"
[[ -x $executable ]] || fail 'Installation did not provide panmanUtils.'
mkdir -p "$prefix/bin"
{
  printf '#!/usr/bin/env bash\nset -euo pipefail\n'
  printf 'exec %q run --no-capture-output -p %q %q "$@"\n' \
    "$forge/bin/conda" "$env_dir" "$executable"
} > "$prefix/bin/panmanUtils"
chmod +x "$prefix/bin/panmanUtils"
"$forge/bin/conda" list -p "$env_dir" --explicit > "$prefix/explicit-spec.txt"
"$prefix/bin/panmanUtils" --help
printf '\nInstalled. Run:\n  %q --help\n' "$prefix/bin/panmanUtils"
