#!/usr/bin/env bash
# Install the released Bioconda package, not the local source checkout.
set -euo pipefail

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
usage() {
  cat <<'HELP'
Usage: bash installationMacOS.sh [--prefix DIR] [--version VERSION]

Install released PanMAN using a private Miniforge environment on macOS 11+.
Default prefix: $HOME/panman (must not contain spaces).
No sudo, shell activation, or shell configuration changes are required.
Run PREFIX/bin/panmanUtils after installation.
This does not compile the current checkout; --version selects a Bioconda release.
HELP
}
prefix="$HOME/panman"
version=''
while (($#)); do
  case "$1" in
    --prefix|--version)
      (($# >= 2)) || fail "$1 requires a value"
      [[ -n $2 && $2 != --* ]] || fail "$1 requires a value"
      case "$1" in --prefix) prefix=$2 ;; --version) version=$2 ;; esac
      shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown option: $1" ;;
  esac
done
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
"$forge/bin/conda" "$operation" -y -p "$env_dir" \
  --override-channels --strict-channel-priority -c conda-forge -c bioconda "$spec"
[[ -x $env_dir/bin/panmanUtils ]] || fail 'The package did not provide panmanUtils.'
mkdir -p "$prefix/bin"
{
  printf '#!/usr/bin/env bash\nset -euo pipefail\n'
  printf 'exec %q run --no-capture-output -p %q %q "$@"\n' \
    "$forge/bin/conda" "$env_dir" "$env_dir/bin/panmanUtils"
} > "$prefix/bin/panmanUtils"
chmod +x "$prefix/bin/panmanUtils"
"$forge/bin/conda" list -p "$env_dir" --explicit > "$prefix/explicit-spec.txt"
"$prefix/bin/panmanUtils" --help
printf '\nInstalled. Run:\n  %q --help\n' "$prefix/bin/panmanUtils"
