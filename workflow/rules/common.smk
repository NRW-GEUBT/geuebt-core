import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read().strip()
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


def validate_path(w):
    """
    returns path to subworklfow directory if param is 'default', returns param otherwise
    """
    if config["geuebt-validate_path"] == "default":
        return f"{str(os.path.expanduser('~'))}/.nrw-geuebt/geuebt-core_{version}/geuebt-validate/workflow/Snakefile"
    return config["geuebt-validate_path"] 


def chewie_path(w):
    """
    returns path to subworklfow directory if param is 'default', returns param otherwise
    """
    if config["geuebt-chewie_path"] == "default":
        return f"{str(os.path.expanduser('~'))}/.nrw-geuebt/geuebt-core_{version}/geuebt-chewie/workflow/Snakefile"
    return config["geuebt-chewie_path"] 


def charak_path(w):
    """
    returns path to subworklfow directory if param is 'default', returns param otherwise
    """
    if config["geuebt-charak_path"] == "default":
        return f"{str(os.path.expanduser('~'))}/.nrw-geuebt/geuebt-core_{version}/geuebt-charak/workflow/Snakefile"
    return config["geuebt-charak_path"] 


# Input functions ------------------------------------
def aggregate_over_species(wildcards):
    "Aggregate chewie and charak rule outputs over the species wildcard and returns a dict of file lists"
    checkpoint_output = checkpoints.create_sample_sheets.get(**wildcards).output["dirout"]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{species}.tsv")).species
    wdict= {
        "chewie_qc": expand(
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
        "charak_qc": expand(
            "charak/{species}/staging/qc_status.json", species=ids_map
        ),
    }
    return(wdict)
