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


USERNAME = os.getenv("GEUEBT_API_USERNAME")
PASSWORD = os.getenv("GEUEBT_API_PASSWORD")


def login(url, username, password):
    response = requests.post(
        f"{url}/users/token",
        data={"username": username, "password": password},
        headers={"Content-Type": "application/x-www-form-urlencoded"}
    )
    response.raise_for_status()
    return response.json()["access_token"]


def authenticated_request( method, endpoint, token, **kwargs):
    headers = kwargs.pop("headers", {})
    headers["Authorization"] = f"Bearer {token}"
    return requests.request(method, endpoint, headers=headers, **kwargs)


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
            {"STATUS": "WARNING", "MESSAGES": ["No geuebt-chewie data for this sample"]}
        )
        
        # merge status
        qc_status = [
            valiqc[k]["STATUS"],
            chewie_value["STATUS"],
            charak_value["STATUS"],
        ]
        if "FAIL" in qc_status:
            valiqc[k]["STATUS"] = "FAIL"
        elif "WARNING"in status:
            valiqc[k]["STATUS"] = "WARNING"
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
    if not USERNAME or not PASSWORD:
        raise RuntimeError("Missing API_USERNAME or API_PASSWORD env vars")
    token = login(url, USERNAME, PASSWORD)
    response = authenticated_request("POST", urljoin(url, "runs"), token , json=qcstatus)
    if response.status_code != 200:
        print(json.dumps(response.json(), indent=4), file=sys.stderr)
        raise ValueError("Failed to POST run status. Check the logs for more details.")


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
