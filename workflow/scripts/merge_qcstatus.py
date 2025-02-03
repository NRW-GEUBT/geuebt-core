#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
import pwd
import os
from urllib.parse import urljoin
import requests


def main(vali_status, chewie_status, charak_status, status, ver, workdir_path, url):
    chewieqc, charakqc = {}, {}
    # merge over species
    for filepath in chewie_status:
        with open(filepath, "r") as fi:
            chewieqc.update(json.load(fi))
    for filepath in charak_status:
        with open(filepath, "r") as fi:
            charakqc.update(json.load(fi))
    # merge validation and calling status
    with open(vali_status, "r") as fi:
        valiqc = json.load(fi)

    for k in valiqc.keys():
        chewie_value = chewieqc.get(
            k, 
            {"STATUS": "FAIL", "MESSAGES": ["No geuebt-chewie data for this sample"]}
        )
        charak_value = charakqc.get(
            k,
            {"STATUS": "WARN", "MESSAGES": ["No geuebt-chewie data for this sample"]}
        )
        
        # merge status
        qc_status = [
            valiqc[k]["STATUS"],
            chewie_value["STATUS"],
            charak_value["STATUS"],
        ]
        if "FAIL" in qc_status:
            valiqc[k]["STATUS"] = "FAIL"
        elif "WARN"in status:
            valiqc[k]["STATUS"] = "WARN"
        else:
            pass
        
        # merge messages
        for msg in chewie_value["MESSAGES"], charak_value["MESSAGES"]:
            valiqc[k]["MESSAGES"].extend(msg)

    # add head info and format
    qcstatus = {
        "run_metadata": {
            "name": os.path.basename(workdir_path),
            "geuebt_version": ver,
            "user": pwd.getpwuid(os.getuid()).pw_name
        }
    }
    qcstatus["samples"] = [
        {
            "isolate_id" : k ,
            "STATUS": v["STATUS"],
            "MESSAGES": v["MESSAGES"]
        } for k, v in valiqc.items()
    ]

    # Write as JSON
    with open(status, "w") as fo:
        json.dump(qcstatus, fo, indent=4)

    # post run data
    response = requests.post(urljoin(url, "runs"), json=qcstatus)
    #Should check repsonse


if __name__ == '__main__':
    main(
        snakemake.input['vali_status'],
        snakemake.input['chewie_status'],
        snakemake.input['charak_status'],
        snakemake.output['status'],
        snakemake.params["geuebt_version"],
        snakemake.params["workdir_path"],
        snakemake.params["url"],
    )
