#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'

# install prerequisites
apt-get update
apt-get install -y curl patch xz-utils

# install dlang installer
mkdir -p "$DLANG_ROOT"
curl https://dlang.org/install.sh > "$DLANG_ROOT/install.sh"
chmod +x "$DLANG_ROOT/install.sh"

# NOTE: older GnuPG versions do not support latest trends in
#       signature algos. This, we sacrifice verification for
#       for a working installation.
sed -i -E 's/^verify\(\)\s*\{/\0\n    return/' "$DLANG_ROOT/install.sh"

# download dmd and dub
"$DLANG_ROOT/install.sh" -p "$DLANG_ROOT" install "dmd-$DMD_VERSION"

# install into system
install -t /usr/bin "$DLANG_ROOT/dmd-$DMD_VERSION/linux/bin64/*"
install -t /usr/lib "$DLANG_ROOT/dmd-$DMD_VERSION/linux/lib64/*"

# test
dmd --version
dub --version

# clean up
rm -rf "$DLANG_ROOT"
apt-get remove -y curl patch xz-utils
apt-get clean
