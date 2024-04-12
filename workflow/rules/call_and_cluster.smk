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
        external_main_clusters="db_data/{species}/clusters.tsv",
        external_sub_clusters="db_data/{species}/subclusters.tsv",
        external_profiles="db_data/{species}/profiles.tsv",
        external_timestamps="db_data/{species}/timestamps.tsv",
        external_statistics="db_data/{species}/statistics.tsv",
    output:
        workdir=directory("call_and_cluster/{species}"),
        qc_status="call_and_cluster/{species}/staging/qc_status.json",
        isolate_sheets="call_and_cluster/{species}/staging/isolate_sheets.json",
        clusters="call_and_cluster/{species}/staging/clusters.json",
    params:
        max_threads_per_job=config["max_threads_per_job"],
        chewie_path=f"{config['geuebt-chewie_path']}/workflow/Snakefile",
        conda_prefix=get_conda_prefix,
        # Using functions to get species specific parameters from config
        prodigal=lambda w: config[w.species]["prodigal"],
        cluster_distance=lambda w: config[w.species]["cluster_distance"],
        subcluster_distance=lambda w: config[w.species]["subcluster_distance"],
        cluster_prefix=lambda w: config[w.species]["cluster_prefix"],
        cgmlst_scheme=lambda w: config[w.species]["cgmlst_scheme"],
        organism=lambda w: f"{w.species.replace('_', ' ').replace('spp', 'spp.')}",
        # Workflow needs absolute paths! And using lambda to get wildcard value
        sample_sheet=lambda w: f"{os.getcwd()}/sample_sheets/{w.species}.tsv",
        external_main_clusters=lambda w: f"{os.getcwd()}/db_data/{w.species}/clusters.tsv",
        external_sub_clusters=lambda w: f"{os.getcwd()}/db_data/{w.species}/subclusters.tsv",
        external_profiles=lambda w: f"{os.getcwd()}/db_data/{w.species}/profiles.tsv",
        external_timestamps=lambda w: f"{os.getcwd()}/db_data/{w.species}/timestamps.tsv",
        external_statistics=lambda w: f"{os.getcwd()}/db_data/{w.species}/statistics.tsv",
    message:
        "[Call and cluster] Allele calling and clustering for {wildcards.species}"
    conda:
        "../envs/git_workflow.yaml"
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
                     external_main_clusters={params.external_main_clusters} \
                     external_sub_clusters={params.external_sub_clusters} \
                     external_profiles={params.external_profiles} \
                     external_timestamps={params.external_timestamps} \
                     external_statistics={params.external_statistics} \
                     prodigal={params.prodigal} \
                     cgmlst_scheme={params.cgmlst_scheme} \
                     max_threads_per_job={params.max_threads_per_job} \
                     max_missing_loci=0.05 \
                     distance_method=3 \
                     clustering_method="single" \
                     cluster_distance={params.cluster_distance} \
                     subcluster_distance={params.subcluster_distance} \
                     cluster_prefix={params.cluster_prefix} \
                     organism="{params.organism}"
        """


rule merge_qcstatus:
    input:
        vali_status="validation/staging/validation_status_ids_checked.json",
        chewie_status=lambda w: aggregate_over_species(w)["qc_status"],
    output:
        status="staging/qc_status.json",
    params:
        geuebt_version=version,
        workdir_path=os.getcwd(),
    message:
        "[Call and cluster] Updating QC status"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_qcstatus.log",
    script:
        "../scripts/merge_qcstatus.py"


checkpoint move_and_split_chewie_results:
    input:
        isolate_sheets=lambda w: aggregate_over_species(w)["isolate_sheets"],
        clusters=lambda w: aggregate_over_species(w)["clusters"],
    output:
        # move clusters to root staging directly since they won't be processed further
        isolates=directory("call_and_cluster/staging/isolates_sheets"),
        clusters=directory("staging/clusters"),
    message:
        "[Call and cluster] Gathering and splitting chewie results"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/move_isolate_sheets_chewie.log",
    script:
        "../scripts/move_and_split_chewie_results.py"
