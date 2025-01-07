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


def main(isolate_sheets, dirout, fasta_prefix):
    if not os.path.isdir(dirout):
        os.mkdir(dirout)
    ssheets = {}
    with open(isolate_sheets, "r") as fi:
        samples = json.load(fi)

    for sample in samples:
        organism = sample["organism"]
        if organism not in ssheets.keys():
            ssheets[organism] = pd.DataFrame(columns=["sample", "assembly"])

        fasta_path = os.path.abspath(
            os.path.join(fasta_prefix, sample["fasta_name"])
        )
        ssheets[organism] = pd.concat(
            [
                ssheets[organism],
                pd.DataFrame({"sample": [sample['isolate_id']], "assembly": [fasta_path]})
            ],
            ignore_index=True
        )
    for k, v in ssheets.items():
        filename = k.replace(" ", "_").replace("spp.", "spp") + ".tsv"
        # v.to_csv(os.path.join(dirout, f"{organism}.tsv"), sep="\t", index=False)
        v.to_csv(os.path.join(dirout, filename), sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input['isolate_sheets'],
        snakemake.output['dirout'],
        snakemake.params['fasta_prefix'],
    )
