# extracting-new-celltypes

These scripts are for testing/extracting cell counts using the new whole blood reference with 12 cell types: [paper-link](https://www.nature.com/articles/s41467-021-27864-7). The code to run this is pulled from the meffil package script [idoloptimized-references.r](https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r).


## Scripts

* [test.R](test.R) is code used to test what format the data needs to be in to run the code
* [extract-cellcounts.R](extract-cellcounts.R) extracts the cell counts from the 15up individuals
* [extract-cellcounts.sh](extract-cellcounts.sh) submits the R script as a job to BC4

Run [extract-cellcounts.sh](extract-cellcounts.sh) using: 

``` bash
sbatch extract-cellcounts.sh
```
