#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import json
import glob
from datetime import datetime


FASTA_EXT = [".fa", ".faa", ".fna", ".ffn", ".frn", ".fasta"]
TABLE_EXT = [".csv", ".tsv", ".txt"]


def iso_timestamp(file):
    return datetime.fromtimestamp(os.path.getmtime(file)).isoformat()


def find_files(folder, exts):
    files = []
    for ext in exts:
        pathstring = os.path.join(folder, f"*{ext}")
        files.extend(glob.glob(pathstring))
    return files


def main(watched_folders, manifest_path, fasta_ext=FASTA_EXT, table_ext=TABLE_EXT):
    # In each watched folder look for fasta files
    filelist = []
    for fd in watched_folders:
        filelist.extend(find_files(fd, fasta_ext))
    fastas = [
        {"path": filepath, "timestamp": iso_timestamp(filepath)}
        for filepath in filelist
    ]
    # look for tables, more than one is allowed
    filelist = []
    for fd in watched_folders:
        filelist.extend(find_files(fd, table_ext))
    tables = [
        {"path": filepath, "timestamp": iso_timestamp(filepath)}
        for filepath in filelist
    ]
    # make manifest as json
    jsonout = {"fasta": fastas, "metadata": tables}
    with open(manifest_path, "w") as fi:
        json.dump(jsonout, fi, indent=4)


if __name__ == '__main__':
    main(
        snakemake.params['watched_folders'],
        snakemake.output['manifest']
    )
