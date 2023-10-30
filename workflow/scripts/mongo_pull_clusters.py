#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import glob
from json import dumps
from pymongo import MongoClient
import pandas as pd


FASTA_EXT = [".fa", ".faa", ".fna", ".ffn", ".frn", ".fasta"]


def db_connect(host, port, database):
    client = MongoClient(host, port)
    db = client[database]
    return db


def list_scheme_files(folder, exts=FASTA_EXT):
    files = []
    for ext in exts:
        pathstring = os.path.join(folder, f"*{ext}")
        files.extend(glob.glob(pathstring))
    return files


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

    for d, header, path in zip(
        [main, sub, timestamps, statistics],
        [
            "sample\tcluster_name\n",
            "sample\tcluster_name\n",
            "sample\tdate\n",
            "sample\tEXC\tINF\tLNF\tPLOT\tNIPH\tALM\tASM\n"
        ],
        [clusters_path, subclusters_path, timestamps_path, statistics_path]
    ):
        # if no data, write the header to file
        if not d:
            with open(path, "w") as fo:
                fo.write(header)
        else:
            # Dump everything in tables
            tbl = pd.read_json(dumps(d), orient="records")
            tbl.to_csv(path, sep="\t", index=False)

    # For the profile, same but need to read in the file names from scheme
    if not profiles:
        tab = "\t"  # backlash forbidden in f-strings
        fastanames = [os.path.split(p)[1] for p in list_scheme_files(scheme)]
        header = f"#FILE\t{tab.join(fastanames)}\n"
        with open(profiles_path, "w") as fo:
                fo.write(header)
    else:
        tbl = pd.read_json(dumps(profiles), orient="records")
        tbl.to_csv(profiles_path, sep="\t", index=False)


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
