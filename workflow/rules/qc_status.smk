#mergeQC status from different workflows and push to DB

rule merge_qcstatus:
    input:
        vali_status="validation/staging/validation_status.json",
        chewie_status=lambda w: aggregate_over_species(w)["chewie_qc"],
        charak_status=lambda w: aggregate_over_species(w)["charak_qc"],
    output:
        status="staging/qc_status.json",
    params:
        geuebt_version=version,
        workdir_path=config["workdir"],
        url=config["API_url"],
    message:
        "[Call and cluster] Updating QC status"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_qcstatus.log",
    script:
        "../scripts/merge_qcstatus.py"
