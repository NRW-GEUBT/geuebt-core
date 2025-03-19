### 1.2.5

Bump chewie to 1.3.2

### 1.2.4

Bump chewie to 1.3.1

### 1.2.3

Correctly cleanup metadata from input dirs

### 1.2.0

Fix date parsing (validate)
Add workflow version to isolate_sheet (all)

### 1.1.1

Adding missing test data

### 1.1.0

Deploy all dependencies under CONDA_PREFIX

### 0.2.10

Modified the deloyement of subworkflows to allow or concurrent versions of core and/or the subworkflows to be used on the same system without clashes
Update validate to 0.2.7

### 0.2.9

Fix grapetree crash with null distance matrix

### 0.2.8

Fix a issue were cluster would not appear in the trees for orphan "clusters"

### 0.2.7

Add trees in Newick format to the cluster sheets
Fix test DB in test configuration

### 0.2.6

Update validate to fix date parsing and better handling of isolates with lots of plasmids

### 0.2.5

Update BakCharak
Known Issue: Workflow is not working with snakemake8 for the time being. Make sure you use snakemake7.

### 0.2.4

Fix funky run order

### 0.2.3

Now register run name (Workdir name) in the run sheet
Now waits for all critical steps to finish before commiting to database
Bump validate to 0.2.5
Bump chewie to 0.3.2

### 0.2.2

Bump geuebt-chewie to 0.3.1

### 0.2.1

Bump geuebt-validate to 0.2.1

### 0.1.8

Bump validate module to 0.1.1

### 0.1.7

- Fix organism name in cluster sheets
- Validate bump to 0.1.7

### 0.1.6

- Fix sample sheet creation when mutliple organisms present in dataset
- Fix creation of empty cluster and isolates desription files when the database is initially empty

### 0.1.5

Bump geuebt-validate to v0.1.6

### 0.1.4

Bump geuebt-validate to v0.1.5

### 0.1.3

Fix generation of files for chewie when the database is initally empty

### 0.1.2

Bump geuebt-chewie to v0.3.0

### 0.1.1

Will reject samples whose field 'isolate_id' is already used in the database.
THis effectivaly happens after qc validation so it will be visible only for those
samples that pass qc.

### 0.1.0

First release
