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


def main(manifest, metadata_path, cleanup):
    with open(manifest, "r") as fi:
        manifest_dict = json.load(fi)
    metadata_imputs = [entry["path"] for entry in manifest_dict["metadata"]]
    dfs = [
        pd.read_csv(path, sep="\t")
        for path in metadata_imputs
    ]
    if len(dfs) > 1:
        merged = pd.concat(dfs, axis=0)
    else:
        merged = dfs[0]
    merged.to_csv(metadata_path, sep="\t", index=False)
    if cleanup:
        for path in metadata_imputs:
            os.remove(path)


if __name__ == '__main__':
    main(
        snakemake.input['manifest'],
        snakemake.output['metadata'],
        snakemake.params['cleanup'],
    )
