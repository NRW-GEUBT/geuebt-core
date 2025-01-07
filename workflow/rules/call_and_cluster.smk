# Split sample by species and perform allele calling and clustering speices-wise
# Samrt way would be to check which species are there and only run these
# rules that need to be run
# HINT: make a template rule wrapping chewie module then make subroules with
# parameters adapted to each species. Read in the different species with a function
# and run only the nescessary rules. At the end aggregate everything to staging


import os


rule chewie:
    # Inputs need absolute paths!
    input:
        sample_sheet="sample_sheets/{species}.tsv",
    output:
        workdir=directory("call_and_cluster/{species}"),
        qc_status="call_and_cluster/{species}/staging/qc_status.json",
        isolate_sheets="call_and_cluster/{species}/staging/isolate_sheets.json",
        clusters="call_and_cluster/{species}/staging/clusters.json",
    params:
        max_threads_per_job=config["max_threads_per_job"],
        chewie_path=os.path.expanduser(f"~/.nrw-geuebt/geuebt-core-{version}/geuebt-chewie/workflow/Snakefile"),
        conda_prefix=get_conda_prefix,
        organism=lambda w: f"{w.species.replace('_', ' ').replace('spp', 'spp.')}",
        url=config["API_url"],
        ephemeral=config["ephemeral"],
        # Workflow needs absolute paths! And using lambda to get wildcard value
        sample_sheet=lambda w: f"{os.getcwd()}/sample_sheets/{w.species}.tsv",
    message:
        "[Call and cluster] Allele calling and clustering for {wildcards.species}"
    conda:
        "../envs/chewie.yaml"
    threads: workflow.cores
    log:
        "logs/call_and_cluster_{species}.log",
    shell:
        """
        exec 2> {log}
        snakemake -s {params.chewie_path} \
            --use-conda \
            --conda-prefix {params.conda_prefix} \
            --cores {threads} \
            --config workdir={output.workdir} \
                     sample_sheet={params.sample_sheet} \
                     API_url={params.url}\
                     max_threads_per_job={params.max_threads_per_job} \
                     ephemeral={params.ephemeral} \
                     organism="{params.organism}"
        """


checkpoint move_and_split_chewie_results:
    input:
        isolate_sheets=lambda w: aggregate_over_species(w)["isolate_sheets"],
        clusters=lambda w: aggregate_over_species(w)["clusters"],
    output:
        isolates=directory("call_and_cluster/staging/isolates_sheets"),
        clusters=directory("call_and_cluster/staging/clusters/"),
    message:
        "[Call and cluster] Gathering and splitting chewie results"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/move_isolate_sheets_chewie.log",
    script:
        "../scripts/move_and_split_chewie_results.py"
