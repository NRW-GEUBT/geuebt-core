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
from pymongo import MongoClient


def db_connect(host, port, database):
    client = MongoClient(host, port)
    db = client[database]
    return db


def main(isolate_sheets, qc_in, dirout, qc_out, fasta_prefix, host, port, database):
    if not os.path.isdir(dirout):
        os.mkdir(dirout)
    with open(qc_in, "r") as fi:
        qc = json.load(fi)
    ssheets = {}
    with open(isolate_sheets, "r") as fi:
        samples = json.load(fi)
    for sample in samples:
        # Checking if isolate_id is already in use
        # if it is then do not add sample to sample sheet and add error message to QC_sheet.
        collection = db_connect(host, port, database)['isolates']
        if collection.find_one({"isolate_id": sample["isolate_id"]}):
            qc[sample["isolate_id"]]["STATUS"] = "FAIL"
            qc[sample["isolate_id"]]["MESSAGES"].append("This isolate_id is already in use.")
            continue

        # if the isolate_id is available then add smaple to ssample sheet
        if not sample["organism"] in ssheets:
            organism = sample["organism"]
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
        filename = k.replace(" ", "_").replace(".", "") + ".tsv"
        v.to_csv(os.path.join(dirout, filename), sep="\t", index=False)
    with open(qc_out, "w") as fo:
        json.dump(qc, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['isolate_sheets'],
        snakemake.input['qc_in'],
        snakemake.output['dirout'],
        snakemake.output['qc_out'],
        snakemake.params['fasta_prefix'],
        snakemake.params['host'],
        snakemake.params['port'],
        snakemake.params['database'],
    )
