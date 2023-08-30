#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import pathlib
from json import dumps
from pymongo import MongoClient
import pandas as pd


def db_connect(host, port, database):
    client = MongoClient(host, port)
    db = client[database]
    return db


def main(
    clusters_path,
    subclusters_path,
    profiles_path,
    timestamps_path,
    statistics_path,
    species_wildcard,
    host,
    port,
    database,
    scheme
):
    species = species_wildcard.replace("_", " ").replace("spp", "spp.")
    db = db_connect(host, port, database)

    # Getting clusters populations
    cluster_coll = db["clusters"]
    main, sub = [], []
    for clus in cluster_coll.find({
        "organism": species},
        comment="Automated <find> of cluster sheets from Geuebt-Core"
    ):
        # Oprhans are skipped
        if clus["cluster_number"] > 0:
            main.append({"sample": clus["representative"], "cluster_name": clus["cluster_id"]})
            for subcluster in clus["subclusters"]:
                sub.append({"sample": subcluster["representative"], "cluster_name": subcluster["subcluster_id"]})

    # Getting sample info
    isolates_coll = db["isolates"]
    profiles, timestamps, statistics = [], [], []
    for sample in isolates_coll.find(
        {"organism": species},
        comment="Automated <find> of cluster sheets from Geuebt-Core"
    ):
        timestamps.append({"sample": sample["isolate_id"], "date": sample["epidata"]["collection_date"]})
        statistics.append({"sample": sample["isolate_id"], **sample["cgmlst"]["allele_stats"]})
        massaged_profiles = {entry["locus"]: entry["allele_crc32"] for entry in sample["cgmlst"]["allele_profile"]}
        profiles.append({"#FILE": sample["isolate_id"], **massaged_profiles})

    # Ensure valid files headers
    for d, header in zip(
        [main, sub, timestamps, statistics],
        [
            {"sample": "", "cluster_name": ""},
            {"sample": "", "cluster_name": ""},
            {"sample": "", "date": ""},
            {"sample": "", "EXC": "", "INF": "", "LNF": "", "PLOT": "", "NIPH": "", "ALM": "", "ASM": ""}
        ]
    ):
        if not d:
            d = header
    # For the profile, need to read in the file names from scheme
    if not profiles:
        fastanames = list(pathlib.Path(scheme).glob('*.fasta'))
        profiles = {"#FILE": ""}
        profiles.update({filename: "" for filename in fastanames})

    # Dump everything in tables
    for d, path in zip(
        [main, sub, timestamps, statistics, profiles],
        [clusters_path, subclusters_path, timestamps_path, statistics_path, profiles_path]
    ):
        tbl = pd.read_json(dumps(d), orient="records")
        tbl.to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    main(
        snakemake.output["clusters"],
        snakemake.output["subclusters"],
        snakemake.output["profiles"],
        snakemake.output["timestamps"],
        snakemake.output["statistics"],
        snakemake.params["organism"],
        snakemake.params["host"],
        snakemake.params["port"],
        snakemake.params["database"],
        snakemake.params["cgmlst_scheme"]
    )
