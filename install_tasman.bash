#!/usr/bin/env bash

set -euo pipefail

# Determine the TASMAN installation directory, independently of where
# the script is called from.

tasman_dir=$(cd "$(dirname "$0")" && pwd)
source_dir="$tasman_dir/source"

# Verify that the expected directories and build files exist.

if [[ ! -d "$source_dir" ]]; then
  echo "TASMAN installation error: source directory not found:" >&2
  echo "  $source_dir" >&2
  exit 1
fi

if [[ ! -f "$source_dir/Makefile" ]]; then
  echo "TASMAN installation error: Makefile not found:" >&2
  echo "  $source_dir/Makefile" >&2
  exit 1
fi

misc_file="$tasman_dir/misc/score.tab"

if [[ ! -f "$misc_file" ]]; then
  echo "TASMAN installation error: misc database missing or incomplete:" >&2
  echo "  $misc_file" >&2
  exit 1
fi

echo
echo "Installing TASMAN"
echo "Installation directory: $tasman_dir"
echo

# Pass all command-line arguments directly to make. This permits, e.g.:
#
# ./install_tasman.bash FC=ifx FFLAGS="-O3"
# ./install_tasman.bash FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"

make -C "$source_dir" clean
make -C "$source_dir" all "$@"

tasman_exe="$tasman_dir/bin/tasman"

if [[ ! -x "$tasman_exe" ]]; then
  echo "TASMAN installation error: executable not created:" >&2
  echo "  $tasman_exe" >&2
  exit 1
fi

echo
echo "TASMAN executable:"
echo "  $tasman_exe"
echo
echo "If not already done, add the following lines to your shell configuration:"
echo
echo "  export TASMAN_DIR=\"$tasman_dir\""
echo "  export PATH=\"\$TASMAN_DIR/bin:\$PATH\""
echo "  export TASMAN_USER=\"Your Name\""
echo
echo "Alternatively, edit code_dir in source/machine.f90 and rebuild TASMAN."
echo
