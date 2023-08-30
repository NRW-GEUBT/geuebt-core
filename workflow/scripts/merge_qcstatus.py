#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(vali_status, chewie_status, status):
    callqc = {}
    # merge over species
    for filepath in chewie_status:
        with open(filepath, "r") as fi:
            callqc.update(json.load(fi))
    # merge validation and calling status
    with open(vali_status, "r") as fi:
        valiqc = json.load(fi)
    for k in valiqc.keys():
        call_values = callqc.get(
            k,
            {"STATUS": "FAIL", "MESSAGES": ["Allele calling not performed"]}
        )
        valiqc[k]["STATUS"] = call_values["STATUS"]
        valiqc[k]["MESSAGES"].extend(call_values["MESSAGES"])
    with open(status, "w") as fo:
        json.dump(valiqc, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['vali_status'],
        snakemake.input['chewie_status'],
        snakemake.output['status']
    )
