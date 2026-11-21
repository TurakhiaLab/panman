#!/usr/bin/env bash
# Build the local source checkout by default; --package installs a Bioconda release.
set -euo pipefail

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
usage() {
  cat <<'HELP'
Usage: bash installationMacOS.sh [--prefix DIR] [--repo DIR] [--threads N]
                               [--package [--version VERSION]]

Build PanMAN from source using private Miniforge dependencies on macOS 11+.
Requires Apple Command Line Tools: xcode-select --install
Default prefix: REPO/build (must not contain spaces).
Source builds go in REPO/build; --prefix controls dependencies and the wrapper.
No sudo, shell activation, or shell configuration changes are required.
Run PREFIX/bin/panmanUtils after installation.
The source checkout defaults to the parent of this script's install directory.
--package installs a released Bioconda binary instead; --version selects its release.
HELP
}
prefix=''
version=''
mode=source
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
threads=2
while (($#)); do
  case "$1" in
    --prefix|--version|--repo|--threads)
      (($# >= 2)) || fail "$1 requires a value"
      [[ -n $2 && $2 != --* ]] || fail "$1 requires a value"
      case "$1" in
        --prefix) prefix=$2 ;; --version) version=$2 ;;
        --repo) repo=$2 ;; --threads) threads=$2 ;;
      esac
      shift 2 ;;
    --package) mode=package; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown option: $1" ;;
  esac
done
[[ $threads =~ ^[1-9][0-9]*$ ]] || fail '--threads must be a positive integer.'
[[ $mode == package || -z $version ]] || fail '--version requires --package.'
if [[ $mode == source ]]; then
  [[ -f $repo/CMakeLists.txt && -f $repo/src/panmanUtils.cpp ]] || fail "Not a PanMAN source checkout: $repo"
  repo=$(cd "$repo" && pwd -P)
  /usr/bin/xcrun --find clang++ >/dev/null 2>&1 || fail 'Install Apple Command Line Tools: xcode-select --install'
fi
prefix=${prefix:-"$repo/build"}
[[ $(uname -s) == Darwin ]] || fail 'This installer requires macOS.'
arch=$(uname -m)
# Prefer native Apple Silicon packages even from a Rosetta terminal.
if [[ $(/usr/sbin/sysctl -n hw.optional.arm64 2>/dev/null || true) == 1 ]]; then
  arch=arm64
fi
case "$arch" in arm64) platform=osx-arm64 ;; x86_64) platform=osx-64 ;; *) fail "Unsupported architecture: $arch" ;; esac
major=$(/usr/bin/sw_vers -productVersion)
[[ ${major%%.*} -ge 11 ]] || fail 'Miniforge requires macOS 11 or newer.'
[[ -z $version || $version =~ ^[0-9][0-9A-Za-z._-]*$ ]] || fail 'Invalid version.'
command -v curl >/dev/null || fail 'curl is required.'
[[ $prefix != *[[:space:]]* ]] || fail 'Choose a prefix without whitespace (Miniforge requirement).'
mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd -P)
[[ $prefix != *[[:space:]]* ]] || fail 'The resolved prefix contains whitespace.'
forge="$prefix/miniforge"
env_dir="$prefix/env"
[[ $mode != source ]] || env_dir="$prefix/source-env"
# Ignore an active Conda environment or a platform override in the calling shell.
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SUBDIR PYTHONHOME PYTHONPATH
export CONDA_SUBDIR="$platform"
if [[ ! -x "$forge/bin/conda" ]]; then
  [[ ! -e $forge ]] || fail "Incomplete Miniforge installation at $forge; use a new prefix."
  download_dir=$(mktemp -d "${TMPDIR:-/tmp}/panman-install.XXXXXX")
  trap 'rm -rf "$download_dir"' EXIT
  # Checksum assets exist only for versioned installer filenames.
  miniforge_version=26.9.4-2
  installer="Miniforge3-${miniforge_version}-MacOSX-$arch.sh"
  url="https://github.com/conda-forge/miniforge/releases/download/$miniforge_version"
  curl -fL --retry 3 "$url/$installer" -o "$download_dir/$installer"
  curl -fL --retry 3 "$url/$installer.sha256" -o "$download_dir/checksum"
  expected=$(awk 'NR == 1 {print $1}' "$download_dir/checksum")
  [[ $expected =~ ^[a-fA-F0-9]{64}$ ]] || fail 'Invalid Miniforge checksum.'
  (cd "$download_dir" && printf '%s  %s\n' "$expected" "$installer" | /usr/bin/shasum -a 256 -c -)
  bash "$download_dir/$installer" -b -p "$forge"
fi
spec=panman
[[ -z $version ]] || spec="panman=$version"
operation=create
if [[ -f $env_dir/conda-meta/history ]]; then
  operation=install
elif [[ -e $env_dir ]]; then
  fail "Not a Conda environment: $env_dir; use a new prefix."
fi
if [[ $mode == source ]]; then
  # PanMAN still uses task_scheduler_init and other pre-oneTBB APIs.
  "$forge/bin/conda" "$operation" -y -p "$env_dir" \
    --override-channels --strict-channel-priority -c conda-forge \
    'cmake>=3.20,<4' make pkg-config 'libboost-devel>=1.84,<1.86' \
    capnproto=1.0.2 'protobuf>=3.20,<3.21' jsoncpp \
    tbb=2020.2 tbb-devel=2020.2 zlib xz
  build="$repo/build"
  "$forge/bin/conda" run --no-capture-output -p "$env_dir" \
    cmake -S "$repo" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
    -DCMAKE_OSX_ARCHITECTURES="$arch" -DCMAKE_PREFIX_PATH="$env_dir" \
    -DCAPNP_EXECUTABLE="$env_dir/bin/capnp" \
    -DCAPNPC_CXX_EXECUTABLE="$env_dir/bin/capnpc-c++"
  "$forge/bin/conda" run --no-capture-output -p "$env_dir" \
    cmake --build "$build" --target panmanUtils --parallel "$threads"
  executable="$build/panmanUtils"
else
  "$forge/bin/conda" "$operation" -y -p "$env_dir" \
    --override-channels --strict-channel-priority -c conda-forge -c bioconda "$spec"
  executable="$env_dir/bin/panmanUtils"
fi
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
