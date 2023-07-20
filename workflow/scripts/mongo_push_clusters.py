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


def main(cluster_paths, host, port, database):
    collection = db_connect(host, port, database)['clusters']
    # Load clusters
    for fp in cluster_paths:
        with open(fp, "r") as fi:
            cluster = json.load(fi)
        collection.replace_one(
            {"cluster_id": cluster["cluster_id"]},
            cluster,
            upsert = True,
            comment="Automated <replace_one> of cluster sheet from Geuebt-Core"
        )


if __name__ == "__main__":
    main(
        snakemake.input["clusters"],
        snakemake.params["host"],
        snakemake.params["port"],
        snakemake.params["database"],
    )
