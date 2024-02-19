# Split samples by species and perfomr characterization


rule charak:
    # Inputs need absolute paths!
    input:
        sample_sheet="sample_sheets/{species}.tsv",
    output:
        workdir=directory("charak/{species}"),
        isolate_sheets="charak/{species}/staging/merged_sheets.json",
    params:
        max_threads_per_job=config["max_threads_per_job"],
        charak_path=f"{config['geuebt-charak_path']}/workflow/Snakefile",
        conda_prefix={workflow.conda_prefix},
        species_tag=lambda w: config[w.species]["charak_tag"],
        sample_sheet=lambda w: f"{os.getcwd()}/sample_sheets/{w.species}.tsv",
    message:
        "[Cahrak] Characterizing samples for {wildcards.species}"
    conda:
        "../envs/git_workflow.yaml"
    threads: workflow.cores
    log:
        "logs/charak_{species}.log",
    shell:
        """
        exec 2> {log}
        snakemake -s {params.charak_path} \
            --use-conda --keep-incomplete \
            --conda-prefix {params.conda_prefix} \
            --cores {threads} \
            --config workdir={output.workdir} \
                     sample_sheet={params.sample_sheet} \
                     max_threads_per_job={params.max_threads_per_job} \
                     species={params.species_tag}
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
