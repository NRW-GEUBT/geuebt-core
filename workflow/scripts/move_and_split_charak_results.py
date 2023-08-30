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


def load_split_dump(json_in, dirout, id_field):
    with open(json_in, "r") as fi:
        merged = json.load(fi)
    for record in merged:
        with open(os.path.join(dirout, f"{record[id_field]}.json"), "w") as fo:
            json.dump(record, fo, indent=4)


def main(isolate_in, isolate_dir):
    os.makedirs(isolate_dir, exist_ok=True)

    for json_file in isolate_in:
        load_split_dump(json_file, isolate_dir, "isolate_id")


if __name__ == '__main__':
    main(
        snakemake.input['isolate_sheets'],
        snakemake.output['isolates'],
    )
