#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


from pymongo import MongoClient
from pathlib import Path


def db_connect(host, port, database):
    client = MongoClient(host, port)
    db = client[database]
    return db


def main(seqfiles, host, port, database):
    collection = db_connect(host, port, database)['sequences']

    docs = []
    for seqf in seqfiles:
        with open(seqf, "r") as fi:
            docs.append({
                "isolate_id": Path(seqf).stem,
                "sequence_type": "fasta",
                "sequence": fi.read()
            })
    collection.insert_many(docs, comment="Automated <insert_many> of fastas from Geuebt-Core")


if __name__ == "__main__":
    main(
        snakemake.input["fastas"],
        snakemake.params["host"],
        snakemake.params["port"],
        snakemake.params["database"],
    )
