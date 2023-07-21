# Guide for users

## Installation

### Conda

Install conda from any distribution, i.e. miniconda.
You can follow the setup guide from the [Bioconda team](https://bioconda.github.io/).

We advise installing the mamaba solver in the base environement to speed up
environments creation.

```bash
conda install mamba -n base -c conda-forge
```

### Run environment

The running environement simply required a recent python version (>= 3.9) and snakemake.
If you followed the steps above just run:

```bash
mamba create -n snakemake snakemake
```

### MongoDB

Install the correct MongoDB distribution for your system by following the 
[installation guide](https://www.mongodb.com/docs/manual/installation/#mongodb-installation-tutorials)

Check the documentaiton for database optimization, management, and security.

### Install module and databases

Download the [latest realease](https://github.com/NRW-GEUBT/geuebt-core/releases/latest)
and unpack it.

If you're feeling brave, clone the repository form Github:

```bash
git clone https://github.com/NRW-GEUBT/geuebt-core
```

All software and databases dependencies, including the validation and chewie modules
will be installed during the first run.
It may take a little time so be patient!

### cgMLST schemes

cgMLST schemes need to be provided by the user.
A good place to start for these is [pubMLST](https://pubmlst.org/).
In general it is recommended to look for schemes that are globally make a
consensus for the species of interest.

Importantly, the scheme must be formatted to be used by [CHEWBBACCA](https://github.com/B-UMMI/chewBBACA).

## Configuration

The configuaration can be defined in two ways:

- either edit and locally save the `config/config.yaml` files and provide its path
  to the snakemake command with the `--configfile` argument

- or provide the parameters directly to the snakemake command with
  `--config <ARGNAME>=<VALUE>`

### User defined parameters

Following arguments must be provided for each run:

| Parameter | Type | Description |
| --- | --- | --- |
| `workdir` | path-like string | Path to the ouptut directory |
| `watched_folders` | list of path-like string | Paths to the folders containing input data (metadata and assemblies) |
| `fasta_store` | path-like string | Paths to the stroage folder for processed assemblies on the file system |

### Organism specific parameters

Organism specific parameters must be provided as a dictionnary whose name is the 
organism name as provided in the metadata with spaces replaced by underscore and trailing 
dot stripped (e.g. `Campylobacter spp.` becomes `Campylobacter_spp`)

Following named parameters are required:

 Parameter | Type | Description |
| --- | --- | --- |
| `prodigal` | path-like string | Path to the prodigal training file for this organism. <br> prodigal training files for many species are >br> provided by the local copy of CHEWBBACCA under <br> "~/.nrw-geuebt/chewieSnake/chewBBACA/CHEWBBACA/prodigal_training_files/" |
| `cgmlst_scheme` | path-like string | Path to the cgMLST scheme |
| `cluster_distance` | integer | Allele difference threshold for main cluster definition |
| `subcluster_distance` | integer | Allele difference threshold for sub-cluster definition |
| `cluster_prefix` | string | A short string to use as prefix for cluster naming |

### Optional parameters

Following parameters are optional and will revert to defaults if not set:

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `max_threads_per_job` | integer | 1 | Max number of threads assigned to a single job |
| `geuebt-validate_path` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the validate module folder |
| `geuebt-chewie_path` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the chewie module folder |
| `cleanup_watched_dirs` | boolean | False | Whether to delete the data from the watched folder |
| `mongodb_host` | string | localhost | MongoDB host |
| `mongodb_port` | integer | 27017 | MongoDB port |
| `mongodb_database` | string | geuebt-test | Name of the database to store and collect isolates and clusters from |
| `validation_min_contig_length` | integer | 500 | Minimal contig length (in bp)<br>to be included in the calculation of assembly<br>statistics |
| `max_missing_loci` | float | 0.05 | Maximumj fraction of loci missing for sample acceptance |
| `distance_method` | int | 3 | Distance calculation method for grapetree |
| `clustering_method` | string | "single" | The agglomeration method to be used for hierarchical <br> clustering. This should be (an unambiguous <br> abbreviation of) one of "ward.D", "ward.D2", "single", <br> "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), <br> "median" (= WPGMC) or "centroid" (= UPGMC) |

## Usage

The workflow can be started with:

```bash
snakemake --use-conda --conda-prefix <PATH TO CONDA ENVS> --configfile <PATH TO CONFIG> --cores <N>
```

See the [snakemake documentaiton](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more information.

## Input formats

### Metadata format

The metadata must be provided as a tab-delimited text data with either of `.txt`, `.csv` or `.tsv` extension,
with strict requirements:

| Field name             | Required                            | Value                                                                                                | Role                                                                                                                                                                            |
| ----------------------- | ------------------------------------ | ----------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| isolate_id             | required                            | Free text (unique)                                                                                   | Unique data identification in database                                                                                                                                          |
| sample_id              | required                            | Free text                                                                                            | Unique sample identification for epidemiological analysis                                                                                                                       |
| organism               | required                            | Choice of ['Salmonella enterica', 'Listeria monocytogenes', 'Escherichia coli', 'Campylobacter spp.'] | Sample categorization for downstream analysis                                                                                                                                   |
| isolate_name_alt       | optional                            | Free text                                                                                            | Laboratory internal name                                                                                                                                                        |
| isolation_org          | required                            | Choice of ['MEL', 'OWL', 'RLD', 'RRW', 'WFL, 'unknown']                                              | Establishing chain of custody and for providing contact information                                                                                                             |
| sequencing_org         | required                            | Choice of ['MEL', 'OWL', 'RLD', 'RRW', 'WFL, 'unknown']                                              | Establishing chain of custody and for providing contact information                                                                                                             |
| extraction_method      | optional                            | Free text                                                                                            | Facilitates the comparison of<br>methodologies, as well as<br>analyses.                                                                                                         |
| library_method         | optional                            | Free text                                                                                            | Facilitates the comparison of<br>methodologies, as well as<br>analyses.                                                                                                         |
| sequencing_instrument  | optional                            | Free text                                                                                            | Facilitates the comparison of<br>methodologies, as well as<br>analyses.                                                                                                         |
| bioinformatics_org     | required                            | Choice of ['MEL', 'OWL', 'RLD', 'RRW', 'WFL, 'unknown']                                              | Establishing chain of custody and for providing contact information                                                                                                             |
| assembly_method        | optional                            | Free text                                                                                            | Facilitates the comparison of<br>methodologies, as well as<br>analyses.                                                                                                         |
| third_party_flag       | required                            | TRUE or FALSE                                                                                        | Mark samples that may be owned by third-party organization                                                                                                                      |
| third_party_owner      | required if third_party_flag is TRUE | Free text                                                                                            | Establishing chain of custody and for providing contact information                                                                                                             |
| sample_type            | required                            | Choice of ['lebensmittel', 'futtermittel', 'tier', 'umfeld', 'human', 'unknown']                     | Sample categorization for downstream analysis                                                                                                                                   |
| fasta_name             | required                            | Free text                                                                                            | Identification of the sequence data associated to this isolate                                                                                                                  |
| fasta_md5              | required                            | Free text                                                                                            | Verification of sequence data integrity                                                                                                                                         |
| collection_date        | required                            | ISO8601 date YYY-MM-DD or 'unknown'                                                                  | Faciliate epidemiological analysis                                                                                                                                              |
| collection_municipality| required                            | Free text (ADV Katalog) or 'unknown'                                                                 | Faciliate epidemiological analysis                                                                                                                                              |
| collection_country     | required                            | ISO3166 alpha-2 code (DE) or 'unknown'                                                               | Faciliate epidemiological analysis                                                                                                                                              |
| collection_cause       | required                            | Free text (ADV Katalog) or 'unknown'                                                                 | Faciliate epidemiological analysis                                                                                                                                              |
| collected_by           | required                            | Free text or 'unknown'                                                                               | Establishing chain of custody and for providing contact information                                                                                                             |
| manufacturer           | optional                            | Free text                                                                                            | Faciliate epidemiological analysis. Can be the manufacturer in case of a food product, sampled production place for environmental samples or animal owner for Veterinary samples |
| designation            | required                            | Free text (ADV Katalog) or 'unknown'                                                                 | Faciliate epidemiological analysis. Description of the type of sample, sampling location or animal host.                                                                        |
| manufacturer_type      | required                            | Free text (ADV Katalog) or 'unknown'                                                                 | Faciliate epidemiological analysis. Type of manufacturer, sampling location or animal husbandery.                                                                               |
| sample_description     | required                            | Free text or 'unknown'                                                                               | Faciliate epidemiological analysis. Procise description of the sample.                                                                                                          |
| lot_number             | optional                            | Free text                                                                                            | Faciliate epidemiological analysis.                                                                                                                                             |
| sequencing_depth       | required                            | Float                                                                                                | Assessing quality and providing a<br>measure of confidence in a<br>sequence.                                                                                                    |
| ref_coverage           | required                            | Float                                                                                                | Assessing quality and providing a<br>measure of confidence in a<br>sequence.                                                                                                    |
| q30                    | required                            | Float                                                                                                | Assessing quality and providing a<br>measure of confidence in a<br>sequence.                                                                                                    |

See an example of metadata table in the `config` folder.

### Sequence files

The assemblies must be provided as fasta files.
Wrapped and unwrapped fastas as well as multifastas are allowed.
There are no special requirements for sequence headers.

## Results

Results to be used for the next steps are located in the `staging` folder in the workdir.

### Status report

A JSON report of the QC checks in the form

```json
{
    "2022-0232977-01": {
        "STATUS": "PASS",
        "MESSAGES": []
    },
    "2022-0232977-02": {
        "STATUS": "FAIL",
        "MESSAGES": [
            "Message checksum",
            "Error message assembly"
        ]
    },
    "2022-0232977-03": {
        "STATUS": "FAIL",
        "MESSAGES": ["Error message checksum"]
    },
    "2022-0232977-04": {
        "STATUS": "FAIL",
        "MESSAGES": [
            "Error message metadata 1",
            "Error message metadata 2"
        ]
    }
}
```

### Isolate metadata and cluster results

Isolate and cluster information are summarized in single JSON files which are pushed to the provided database.

### Assemblies

Assemblies are renamed so that the fastas are in the form `<ISOLATE ID>.fa` and stored in the path
provided by the user.
