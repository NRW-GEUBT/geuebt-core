#!/usr/bin/env bash
set -Eeu


deploy () {
  # if already exists, wipe it clean and redo clone
  [ -d "$local_dir" ] && rm -rf "$local_dir"

  # clone and checkout repo
  echo "Cloning repo and checking out stable tag"
  git clone -q "$repo" "$local_dir"
  cd "$local_dir"
  git checkout "$tag"
}


## Geubt Validate
# Repo URL
repo="https://github.com/NRW-GEUBT/geuebt-validate"

# Commit hash to use
tag="0.1.3"

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/geuebt-validate/"

echo "Deploying Geuebt-validate" && deploy


## Geubt Chewie
# Repo URL
repo="https://github.com/NRW-GEUBT/geuebt-chewie"

# Commit hash to use
tag="0.2.0"

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/geuebt-chewie/"

echo "Deploying Geuebt-chewie" && deploy