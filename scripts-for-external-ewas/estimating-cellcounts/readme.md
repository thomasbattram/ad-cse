# Estimating cell counts

These scripts are for estimating cell counts using the new whole blood reference with 12 cell types: [paper-link](https://www.nature.com/articles/s41467-021-27864-7). 

## Before running any scripts please make a "data" and "results" folder

In linux:

``` bash
mkdir data
mkdir results
```

## Scripts

* [estimate-cellcounts.R](scripts/estimate-cellcounts.R) estimates the cell counts for individuals for using IDAT files
* [estimate-cellcounts.sh](scripts/estimate-cellcounts.sh) submits the above R script as a job to a slurm server - this needs to be updated depending on the server type and file paths
* [estimate-cellcounts-epidish.R](scripts/estimate-cellcounts-epidish.R) estimates cell counts using the recommended EpiDISH procedure
* [compare-cellcounts.R](compare-cellcounts.R) Generates plots to compare the cell counts estimated using meffil and using EpiDISH

To run [estimate-cellcounts.sh](scripts/extract-cellcounts.sh), first you need to edit the script on BC4 and add in the paths missing. Then run the code below:

``` bash
sbatch estimate-cellcounts.sh
```
