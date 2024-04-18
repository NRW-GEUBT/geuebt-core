#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
from datetime import datetime
import pwd
import os


def main(vali_status, chewie_status, status, ver, workdir_path):
    callqc = {}
    # merge over species
    for filepath in chewie_status:
        with open(filepath, "r") as fi:
            callqc.update(json.load(fi))
    # merge validation and calling status
    with open(vali_status, "r") as fi:
        tmp_dict = json.load(fi)
    for k in tmp_dict.keys():
        call_values = callqc.get(
            k,
            {"STATUS": "FAIL", "MESSAGES": ["Allele calling not performed"]}
        )
        tmp_dict[k]["STATUS"] = call_values["STATUS"]
        tmp_dict[k]["MESSAGES"].extend(call_values["MESSAGES"])

    # rearrange entrys in DB friendly manner
    qcstatus = {
        "run_metadata": {
            "name": os.path.basename(workdir_path),
            "date": datetime.now().isoformat(),
            "geuebt_version": ver,
            "user": pwd.getpwuid(os.getuid()).pw_name
        }
    }
    qcstatus["samples"] = [
        {
            "isolate_id" : k ,
            "STATUS": v["STATUS"],
            "MESSAGES": v["MESSAGES"]
        } for k, v in tmp_dict.items()
    ]

    with open(status, "w") as fo:
        json.dump(qcstatus, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['vali_status'],
        snakemake.input['chewie_status'],
        snakemake.output['status'],
        snakemake.params["geuebt_version"],
        snakemake.params["workdir_path"],
    )
