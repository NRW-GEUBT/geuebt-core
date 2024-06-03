#!/usr/bin/env bash
set -Eeu

# Tags
validate_ver="0.2.6"
chewie_ver="0.3.2"
charak_ver="0.2.0"

HERE=$(dirname "$0")
VERSION=$(cat "../../VERSION")

namever_str="geuebt-core_$VERSION"
dirpath="${HOME}/.nrw-geuebt/${namever_str}"

deploy () {
  local repo="$1"
  local tag="$2"
  # Download release
  echo "Downloading release ${tag}"
  wget "https://github.com/NRW-GEUBT/${repo}/archive/refs/tags/${tag}.tar.gz"
  echo "Extracting"
  tar -xzf "${tag}.tar.gz" && mv "${repo}-${tag}" "${repo}"
}

# If directory already exists and this scripts executes then it needs to be refreshed
# Wiping it clean and refetching releases
echo "Deploying subworkflows to ${dirpath}"
[ -d "${dirpath}" ] && echo "${dirpath} already exists, removing..." && rm -rf "${dirpath}"
mkdir -p ${dirpath} && cd ${dirpath}

echo "Deploying geuebt-validate"
deploy "geuebt-validate" "${validate_ver}"

echo "Deploying geuebt-chewie"
deploy "geuebt-chewie" "${chewie_ver}"

echo "Deploying geuebt-charak"
deploy "geuebt-charak" "${charak_ver}"

echo "Cleaning up archives"
rm ./*.tar.gz

exit 0