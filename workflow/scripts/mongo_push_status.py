#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
from pymongo import MongoClient


def db_connect(host, port, database):
    client = MongoClient(host, port)
    db = client[database]
    return db


def main(status, host, port, database):
    collection = db_connect(host, port, database)['runs']
    # Load status
    with open(status, "r") as fi:
        doc = json.load(fi)
    collection.insert_one(doc, comment="Automated <insert_one> of qc_status from Geuebt-Core")


if __name__ == "__main__":
    main(
        snakemake.input["qc_status"],
        snakemake.params["host"],
        snakemake.params["port"],
        snakemake.params["database"],
    )
