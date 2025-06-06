#!/usr/bin/env bash
set -Eeu

# Tags
repover="1.3.3"
reponame="geuebt-charak"

dirpath="${CONDA_PREFIX}/${reponame}"
dirroot=$(dirname $dirpath)

CONDADIR=$(dirname ${CONDA_PREFIX})

# If directory already exists and this scripts executes then it needs to be refreshed
# Wiping it clean and refetching releases
echo "Deploying subworkflows to ${dirpath}"
[ -d "${dirpath}" ] && echo "${dirpath} already exists, removing..." && rm -rf "${dirpath}"
mkdir -p ${dirroot} && cd ${dirroot}

echo "Pulling ${reponame}" >> deploylog.txt
echo "Downloading release ${repover}" >> deploylog.txt
wget "https://github.com/NRW-GEUBT/${reponame}/archive/refs/tags/${repover}.tar.gz" &>> deploylog.txt
echo "Extracting" >> deploylog.txt
tar -xzf "${repover}.tar.gz" &>> deploylog.txt
mv "${reponame}-${repover}" "${reponame}" &>> deploylog.txt
rm "${repover}.tar.gz" &>> deploylog.txt

echo "Deploying ${reponame}" >> deploylog.txt
cd ${dirpath}
# snakemake --use-conda --conda-prefix ${CONDADIR} --cores 1 &>> ${dirroot}/deploylog.txt

exit 0