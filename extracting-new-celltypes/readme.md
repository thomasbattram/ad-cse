# extracting-new-celltypes

These scripts are for testing/extracting cell counts using the new whole blood reference with 12 cell types: [paper-link](https://www.nature.com/articles/s41467-021-27864-7). The code to run this is pulled from the meffil package script [idoloptimized-references.r](https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r).

## To do before running scripts

Before running the scripts follow these steps:

1. use the `aries` R package to identify where the IDAT files are:

``` R
# devtools::install_github("MRCIEU/aries")
library(aries)

aries_dir <- "" ## FILL THIS IN
aries.time.points(aries_dir)

aries <- aries.select(aries_dir, time.point = "F24")
samples <- aries$samples
head(samples)
```

2. Check you have access to IDAT files
3. Check the size of the directory containing the IDAT files using `du -sh PATH_TO_IDAT_FILES`
4. Move IDAT files to a space where they can be used on compute node - i.e. NOT the RDSF

## Scripts

* [test.R](test.R) is code used to test what format the data needs to be in to run the code
* [get-f24-measured-cells.R](get-f24-measured-cells.R) extracts the measured cell numbers from the F24 timepoint
* [extract-cellcounts.R](extract-cellcounts.R) extracts the cell counts for ARIES individuals for a given timepoint using IDAT files
* [extract-cellcounts.sh](extract-cellcounts.sh) submits the above R script as a job to BC4 - this needs to be updated depending on the timepoint you are extracting cell counts for
* [compare-cellcounts.R](compare-cellcounts.R) compares cell counts from current reference and previously well-known reference for the 15up timepoint
* [compare-f24-cellcounts.R](compare-f24-cellcounts.R) compares cell counts from current reference and the measured cell numbers for the F24 timepoint

To run [extract-cellcounts.sh](extract-cellcounts.sh), first you need to edit the script on BC4 and add in the paths missing. Then run the code below:

``` bash
sbatch extract-cellcounts.sh
```
