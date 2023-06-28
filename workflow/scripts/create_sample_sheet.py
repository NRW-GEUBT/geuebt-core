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
        if not sample["organism"] in ssheets:
            ssheets[organism] = pd.DataFrame(columns=["sample","path"])
        fasta_path = os.path.abspath(
            os.path.join(fasta_prefix, sample["fasta_name"]
        )
        ssheets[organism] = ssheets[organism].append(
            {"sample": sample['isolate_id'], "path": fasta_path},
            ignore_index=True,
            inplace=True
        )
    for k, v in ssheets.items():
        v.to_csv(os.path.join(dirout, k), sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input['isolate_sheets'],
        snakemake.output['dirout'],
        snakemake.params['fasta_prefix']
    )
