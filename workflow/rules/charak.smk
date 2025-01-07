# Split samples by species and perfomr characterization


import os


rule charak:
    # Inputs need absolute paths!
    input:
        sample_sheet="sample_sheets/{species}.tsv",
    output:
        workdir=directory("charak/{species}"),
        isolate_sheets="charak/{species}/staging/merged_sheets.json",
        qc="charak/{species}/staging/qc_status.json",
    params:
        max_threads_per_job=config["max_threads_per_job"],
        charak_path=os.path.expanduser(f"~/.nrw-geuebt/geuebt-core-{version}/geuebt-charak/workflow/Snakefile"),
        url=config["API_url"],
        ephemeral=config["ephemeral"],
        conda_prefix=get_conda_prefix,
        species=lambda w: f"{w.species.replace('_', ' ').replace('spp', 'spp.')}",
        sample_sheet=lambda w: f"{os.getcwd()}/sample_sheets/{w.species}.tsv",
    message:
        "[Charak] Characterizing samples for {wildcards.species}"
    conda:
        "../envs/charak.yaml"
    threads: workflow.cores
    log:
        "logs/charak_{species}.log",
    shell:
        """
        exec 2> {log}
        snakemake -s {params.charak_path} \
            --use-conda \
            --conda-prefix {params.conda_prefix} \
            --cores {threads} \
            --config workdir={output.workdir} \
                     sample_sheet={params.sample_sheet} \
                     max_threads_per_job={params.max_threads_per_job} \
                     species="{params.species}" \
                     API_url={params.url} \
                     ephemeral={params.ephemeral}
        """


checkpoint move_and_split_charak_results:
    input:
        isolate_sheets=lambda w: aggregate_over_species(w)["charak_sheets"],
    output:
        # move clusters to root staging directly since they won't be processed further
        isolates=directory("charak/staging/isolates_sheets"),
    message:
        "[Call and cluster] Gathering and splitting charak results"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/move_isolate_sheets_charak.log",
    script:
        "../scripts/move_and_split_charak_results.py"
