# Set of CRUD operation for mongoDB


rule mongo_pull_clusters:
    output:
        clusters="db_data/{species}/clusters.tsv",
        subclusters="db_data/{species}/subclusters.tsv",
        profiles="db_data/{species}/profiles.tsv",
        timestamps="db_data/{species}/timestamps.tsv",
        statistics="db_data/{species}/statistics.tsv",
    params:
        organism=lambda w: w.species,
        host=config["mongodb_host"],
        port=config["mongodb_port"],
        database=config["mongodb_database"],
        cgmlst_scheme=lambda w: config[w.species]["cgmlst_scheme"],
    conda:
        "../envs/mongodb.yaml"
    log:
        "logs/mongo_pull_clusters_{species}.log",
    script:
        "../scripts/mongo_pull_clusters.py"


rule make_isolate_sheet:
    input:
        validation="validation/staging/isolates_sheets/{isolate}.json",
        calling="call_and_cluster/staging/isolates_sheets/{isolate}.json",
        charak="charak/staging/isolates_sheets/{isolate}.json",
    output:
        isolate_sheet="staging/isolates_sheets/{isolate}.json",
    message:
        "[Database operations] Creating isolate sheet"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/make_isolate_sheet_{isolate}.log",
    script:
        "../scripts/make_isolate_sheet.py"


rule merge_qcstatus:
    input:
        vali_status="validation/staging/validation_status_ids_checked.json",
        chewie_status=lambda w: aggregate_over_species(w)["qc_status"],
    output:
        status="staging/qc_status.json",
    params:
        geuebt_version=version,
        workdir_path=config['workdir'],
    message:
        "[Call and cluster] Updating QC status"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_qcstatus.log",
    script:
        "../scripts/merge_qcstatus.py"


rule push_fastas:
    input:
        fastas=aggregate_fastas,
    output:
        flag=touch("dbops/push_fastas.flag"),
    params:
        host=config["mongodb_host"],
        port=config["mongodb_port"],
        database=config["mongodb_database"],
    conda:
        "../envs/mongodb.yaml"
    log:
        "logs/push_fastas.log",
    script:
        "../scripts/mongo_push_fastas.py"


rule mongo_push_clusters:
    input:
        clusters=aggregate_clusters,
    output:
        flag=touch("dbops/push_clusters.flag"),
    params:
        host=config["mongodb_host"],
        port=config["mongodb_port"],
        database=config["mongodb_database"],
    conda:
        "../envs/mongodb.yaml"
    log:
        "logs/mongo_push_clusters.log",
    script:
        "../scripts/mongo_push_clusters.py"


rule mongo_push_isolate:
    input:
        isolates=aggregate_isolate_sheets,
    output:
        flag=touch("dbops/push_isolates.flag"),
    params:
        host=config["mongodb_host"],
        port=config["mongodb_port"],
        database=config["mongodb_database"],
    conda:
        "../envs/mongodb.yaml"
    log:
        "logs/mongo_push_isolates.log",
    script:
        "../scripts/mongo_push_isolates.py"


rule mongo_push_status:
    input:
        qc_status="staging/qc_status.json",
    output:
        flag=touch("dbops/push_status.flag"),
    params:
        host=config["mongodb_host"],
        port=config["mongodb_port"],
        database=config["mongodb_database"],
    conda:
        "../envs/mongodb.yaml"
    log:
        "logs/mongo_push_status.log",
    script:
        "../scripts/mongo_push_status.py"
