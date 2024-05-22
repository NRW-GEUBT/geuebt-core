# Create Trees for each cluster
# look for all samples in cluster
# (for orphan clusters, pull representative sample of each cluster)
# pull allele scheme for each sample ( some are in the pipeline, others are in database)
# use grapetree to generate trees
# append tree to cluster sheet


rule get_all_profiles_in_cluster:
    input:
        cluster = "call_and_cluster/staging/clusters/{cluster}.json",
    output:
        profiles = "trees/{cluster}_profiles.tsv",
    params:
        host=config["mongodb_host"],
        port=config["mongodb_port"],
        database=config["mongodb_database"],
        isolates_dir="call_and_cluster/staging/isolates_sheets",
    conda:
        "../envs/mongodb.yaml"
    message:
        "[Trees] Collecting profiles for {wildcards.cluster}"
    log:
        "logs/{cluster}_get_all_profiles_in_cluster.log",
    script:
        "../scripts/get_all_profiles_in_cluster.py"


rule make_tree:
    input:
        profile = "trees/{cluster}_profiles.tsv",
    output:
        tree = "trees/{cluster}.tre",
    conda:
        "../envs/grapetree.yaml"
    message:
        "[Trees] Generating tree for {wildcards.cluster}"
    log:
        "logs/{cluster}_make_tree.log",
    shell:
        """
        exec 2> {log}
        grapetree -p {input.profile} -m MSTreeV2 > {output.tree}
        """


rule add_tree_to_cluster_sheet:
    input:
        tree = "trees/{cluster}.tre",
        cluster = "call_and_cluster/staging/clusters/{cluster}.json",
    output:
        cluster = "staging/clusters/{cluster}.json",
    conda:
        "../envs/pandas.yaml"
    message:
        "[Trees] Finalizing cluster sheet for {wildcards.cluster}"
    log:
        "logs/{cluster}_add_tree_to_cluster_sheet.log"
    script:
        "../scripts/add_tree_to_cluster_sheet.py"
