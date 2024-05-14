#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(tree_in, json_in, json_out):
    with open(json_in, "r") as fi, open(tree_in, "r") as tree:
        cluster = json.load(fi)
        cluster["tree"] = tree.read()
    with open(json_out,"w") as fo:
        json.dump(cluster, fo, indent=4)
    

if __name__ == '__main__':
    main(
        snakemake.input['tree'],
        snakemake.input['cluster'],
        snakemake.output['cluster']
    )
