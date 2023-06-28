# Split sample by species and perform allele calling and clustering speices-wise
# Samrt way would be to check which species are there and only run these 
# rules that need to be run
# HINT: make a template rule wrapping chewie module then make subroules with 
# parameters adapted to each species. Read in the different species with a function
# and run only the nescessary rules. At the end aggregate everything to staging





checkpoint chewie:
    input:
        sample_sheet="sample_sheets/{species}_samples.tsv",
    output:
        cluster_sheets
        