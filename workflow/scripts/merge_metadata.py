#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
import pandas as pd


def main(manifest, metadata_path):
    with open(manifest, "r") as fi:
        manifest_dict = json.load(fi)
    dfs = [
        pd.read_csv(entry["path"], sep="\t") 
        for entry in manifest_dict["metadata"]
    ]
    merged = pd.concat(dfs, axis=0)
    merged.to_csv(metadata_path, sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input['manifest'],
        snakemake.output['metadata']
    )
