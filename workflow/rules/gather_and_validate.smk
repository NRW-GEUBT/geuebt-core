# Gather inputs for the watched foldersand validate data


import os


rule create_manifest:
    output:
        manifest="tracking/manifest.json",
    params:
        watched_folders=config["watched_folders"],
    message:
        "[Gather and validate] Writting file Manifest"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/create_manifest.log",
    script:
        "../scripts/create_manifest.py"


rule merge_metadata:
    input:
        manifest="tracking/manifest.json",
    output:
        metadata="inputs/metadata.tsv",
    message:
        "[Gather and validate] Gathering metadata"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_metadata.log"
    script:
        "../scripts/merge_metadata.py"


rule gather_fastas:
    input:
        manifest="tracking/manifest.json",
    output:
        fastadir=directory("inputs/fastas"),
        fastaflag=touch("inputs/fasta_flag"),
    params:
        cleanup=config["cleanup_watched_dirs"],
    message:
        "[Gather and validate] Gathering fastas"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/gather_fastas.log",
    script:
        "../scripts/gather_fastas.py"


rule validate_input:
    # Uses the geuebt validate workflow
    # afaik not possible to define module inputs, or too messy
    # so workaround by getting the module form github via env post-deploy script
    # and just calling snakemake
    # Inputs need absolute paths!
    input:
        fastaflag="inputs/fasta_flag",
        metadata="inputs/metadata.tsv",
    output:
        workdir=directory("validation/"),
        isolate_sheets=directory("validation/staging/isolates_sheets"),
        fastas=directory("validation/staging/fastas"),
        qc_status="validation/staging/validation_status.json",
        merged_isolate_sheet="validation/staging/isolates_datasheet.json",
    params:
        max_threads_per_job=config["max_threads_per_job"],
        geva_path=f"{config['geuebt-validate_path']}/workflow/Snakefile",
        conda_prefix={workflow.conda_prefix},
        # Absolute paths needed for workflow
        fastadir=f"{os.getcwd()}/inputs/fastas",
        metadata=f"{os.getcwd()}/inputs/metadata.tsv",
    message:
        "[Gather and validate] Validating user input"
    conda:
        "../envs/git_workflow.yaml"
    threads: 
        workflow.cores
    log:
        "logs/validate_input.log",
    shell:
        """
        exec 2> {log}
        snakemake -s {params.geva_path} \
            --use-conda \
            --conda-prefix {params.conda_prefix} \
            --cores {threads} \
            --config workdir={output.workdir} \
                     metadata={params.metadata} \
                     fasta_dir={params.fastadir} \
                     max_threads_per_job={params.max_threads_per_job} \
                     min_contig_length=500
        """


checkpoint create_sample_sheets:
    input:
        isolate_sheets="validation/staging/isolates_datasheet.json",
    output:
        # files are named "Genus_species.tsv", points are striped
        dirout=directory("sample_sheets"),
    params:
        fasta_prefix="validation/staging/fastas",
    conda:
        "../envs/pandas.yaml"
    message:
        "[Gather and validate] Creating species-wise sample sheets"
    log:
        "logs/create_sample_sheets.log"
    script:
        "../scripts/create_sample_sheet.py"
