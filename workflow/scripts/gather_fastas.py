#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import shutil
import json


def main(manifest, cleanup, fastadir):
    if not os.path.isdir(fastadir):
        os.mkdir(fastadir)
    # if cleaningup then move the files, otherwise copy
    func = shutil.move if cleanup else shutil.copy
    with open(manifest, "r") as fi:
        jdict = json.load(fi)
    for entry in jdict["fasta"]:
        func(
            entry["path"], 
            os.path.join(fastadir, os.path.basename(entry["path"]))
        )


if __name__ == '__main__':
    main(
        snakemake.input['manifest'],
        snakemake.params['cleanup'],
        snakemake.output['fastadir'],
    )
