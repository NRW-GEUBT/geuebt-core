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


def main(cluster_in, profiles_out, host, port, database, isolates_dir, cluster_dir):
    # Get isolates
    with open(cluster_in, "r") as fi:
        cluster = json.load(fi)
    isolates = list(cluster.get("distance_matrix")[0].keys())
    # try:
    #     isolates.extend([member for subcluster in cluster["subclusters"] for member in subcluster["members"]])
    # except KeyError:
    #     pass

    # get local profiles
    local_profiles = []
    for isolate in isolates:
        try:
            with open(os.path.join(isolates_dir, f"{isolate}.json"), "r") as fi:
                d = json.load(fi)
            local_profiles.append({"isolate_id": isolate, "profile": d["cgmlst"]["allele_profile"]})
        except FileNotFoundError:
            pass

    # Get profiles for isolates in db
    db = db_connect(host, port, database)
    profiles = list(db["isolates"].aggregate([
        {"$match": {"isolate_id": {"$in": isolates}}},
        {"$project": {"_id": 0, "isolate_id": 1, "profile": "$cgmlst.allele_profile"}}
    ]))

    # Isolates without profiles should be cluster references, get the representative profile for these
    notfound = set(
        isolates
    ) - set(
        [pro["isolate_id"] for pro in profiles]
    ) - set(
        [pro["isolate_id"] for pro in local_profiles]
    )
    # Add new (local) clusters
    local_clusters = []
    for cluster in notfound:
        try:
            with open(os.path.join(cluster_dir, f"{cluster}.json"), "r") as fi:
                d = json.load(fi)
            local_clusters.append({"cluster_id": cluster, "representative": d["representative"]})
        except FileNotFoundError:
            pass
    for c in local_clusters:
        notfound.remove(c["cluster_id"])

    # Get remaining clusters from DB
    repr_samples = list(db["clusters"].aggregate([
        {"$match": {"cluster_id": {"$in": list(notfound)}}},
        {"$project": {"_id": 0, "cluster_id": 1, "representative": 1}}
    ]))

    repr_mapping = {entry["representative"]: entry["cluster_id"] for entry in repr_samples + local_clusters}
    # Local profiles
    repr_local = []
    for isolate in repr_mapping.keys():
        try:
            with open(os.path.join(isolates_dir, f"{isolate}.json"), "r") as fi:
                d = json.load(fi)
            repr_local.append({"isolate_id": isolate, "profile": d["cgmlst"]["allele_profile"]})
        except FileNotFoundError:
            pass
    # Remote profiles
    repr_profiles = list(db["isolates"].aggregate([
        {"$match": {"isolate_id": {"$in": list(repr_mapping.keys())}}},
        {"$project": {"_id": 0, "isolate_id": 1, "profile": "$cgmlst.allele_profile"}}
    ]))
    all_repr = repr_profiles + repr_local
    for entry in all_repr:
        id = entry["isolate_id"]
        entry["isolate_id"] = repr_mapping[id]

    # format and dump as tsv
    merged_profiles = local_profiles + all_repr + profiles
    reformat = {}
    for entry in merged_profiles:
        reformat.update(
            {entry["isolate_id"]: {locus["locus"]: locus["allele_crc32"] for locus in entry["profile"]}}
        )
    df = pd.read_json(json.dumps(reformat), orient="index")
    df = df.reset_index(names="#FILE")
    df.to_csv(profiles_out, sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input['cluster'],
        snakemake.output['profiles'],
        snakemake.params['host'],
        snakemake.params['port'],
        snakemake.params['database'],
        snakemake.params['isolates_dir'],
        snakemake.params['cluster_dir'],
    )
