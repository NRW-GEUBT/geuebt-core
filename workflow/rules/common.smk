import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config and inputs -----------------------
validate(config, schema="../schema/config.schema.yaml")


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


def get_conda_prefix(wildcards):
    try:
        # snakemake < 8.0
        return workflow.conda_prefix
    except:
        # snakemake > 8
        return workflow.deployment_settings.conda_prefix


# Input functions ------------------------------------
def aggregate_over_species(wildcards):
    "Aggregate chewie and charak rule outputs over the species wildcard and returns a dict of file lists"
    checkpoint_output = checkpoints.create_sample_sheets.get(**wildcards).output[0]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{species}.tsv")).species
    return {
        "qc_status": expand(
            "call_and_cluster/{species}/staging/qc_status.json", species=ids_map
        ),
        "isolate_sheets": expand(
            "call_and_cluster/{species}/staging/isolate_sheets.json", species=ids_map
        ),
        "clusters": expand(
            "call_and_cluster/{species}/staging/clusters.json", species=ids_map
        ),
        "charak_sheets": expand(
            "charak/{species}/staging/merged_sheets.json", species=ids_map
        ),
    }


def aggregate_clusters(wildcards):
    checkpoint_output = checkpoints.move_and_split_chewie_results.get(
        **wildcards
    ).output["clusters"]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{cluster}.json")).cluster
    return expand("staging/clusters/{cluster}.json", cluster=ids_map)


def aggregate_isolate_sheets(wildcards):
    checkpoint_output = checkpoints.move_and_split_chewie_results.get(
        **wildcards
    ).output["isolates"]
    # force claculation of charak checkpoint
    checkpoints.move_and_split_charak_results.get(**wildcards).output["isolates"]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{isolate}.json")).isolate
    return expand("staging/isolates_sheets/{isolate}.json", isolate=ids_map)


def aggregate_fastas(wildcards):
    checkpoint_output = checkpoints.move_and_split_chewie_results.get(
        **wildcards
    ).output["isolates"]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{isolate}.json")).isolate
    return expand("validation/staging/fastas/{isolate}.fa", isolate=ids_map)
