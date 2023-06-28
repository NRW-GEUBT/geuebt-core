#!/usr/bin/env bash
set -Eeu

# Repo URL
validate_repo="https://github.com/NRW-GEUBT/geuebt-validate"

# Commit hash to use
validate_tag="0.1.1"

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/geuebt-validate/"

# if already exists, wipe it clean and redo clone
[ -d "$local_dir" ] && rm -rf "$local_dir"

# clone and checkout repo
echo "Cloning Validate and checking out stable tag"
git clone -q "$validate_repo" "$local_dir"
cd "$local_dir"
git checkout "$validate_tag"