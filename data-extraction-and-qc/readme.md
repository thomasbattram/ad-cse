# Data extraction and QC for ARIES data

Aim: To define AD in ARIES participants at age 15up (and age 7?) and output the data to the HPC

Each script must be run manually here and you need access to ALSPAC data to run these scripts

Scripts:

* [extract-aries-samples.R](scripts/extract-aries-samples.R) : Get ARIES data so that the sample sizes of different AD definitions can be checked amongst the EWAS participants
* [extract-pheno-data.R](scripts/extract-pheno-data.R) : Define AD in ALSPAC, extract covariate data, output the data, and move it to the HPC
* [clean-meth.R](scripts/clean-meth.R) : QC the DNAm data
* []() : Generate PCs
* [clean-pheno.R](scripts/clean-pheno.R) : Clean the phenotype data and only extract the data that will be used in the EWAS
