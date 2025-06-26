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
import pandas as pd


def decode_csv(path):
    """
    Attempt to read_csv with utf-8 encoding and falls back to ANSI if it fails
    Returns a dtaframe
    """
    try:
        return pd.read_csv(path, sep="\t", encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep="\t", encoding="cp1252")


def main(manifest, metadata_path, cleanup):
    with open(manifest, "r") as fi:
        manifest_dict = json.load(fi)
    metadata_inputs = [entry["path"] for entry in manifest_dict["metadata"]]
    dfs = [
        decode_csv(path)
        for path in metadata_inputs
    ]
    if len(dfs) > 1:
        merged = pd.concat(dfs, axis=0)
    else:
        merged = dfs[0]
    merged.to_csv(metadata_path, sep="\t", index=False)
    if cleanup:
        for path in metadata_inputs:
            os.remove(path)


if __name__ == '__main__':
    main(
        snakemake.input['manifest'],
        snakemake.output['metadata'],
        snakemake.params['cleanup'],
    )
