# Identifying cell-type specific associations between DNA methylation and atopic dermatitis - scripts for running analyses outside of Bristol

Scripts are provided to run the following analyses: 

1. Estimation of the 12 blood cell counts 
2. Surrogate variable generation 
3. Conventional EWAS 
4. Cell-type specific EWAS using CellDMC, and TCA 
5. Variance EWAS 
6. Cell-type specific EWAS using omicWAS (to be completed after the meta-analysis has been finalised)

The analyses should be performed by first using your own scripts to QC your DNAm data and phenotype data, as described in the analysis plan sent. Then you should estimate the 12 blood cell counts within your dataset. The code to do this, is within [`estimating-cellcounts`](estimating-cellcounts). Finally, all EWAS (except using omicWAS) should be conducted using the scripts in [`ewas`](ewas). Once the meta-analysis has been completed the summary statistics will be sent back and these can be used to run the EWAS with omicWAS using the scripts in [`omicwas-rep`](omicwas-rep).

Within each of the folders, you'll find another `readme` detailing the workflow in that folder.

The [`ewas`](ewas) analyses can be run in a pipelined manner using the `Snakemake`, can be run manually in R, or can be submitted as individual scripts to be run on a remote server at your institution. Using `Snakemake` should be the easiest way to run the analyses, but you'll need to install `Python`, `conda`, `Snakemake`, and `tmux` to do this. 

* Installing Python: https://python.land/installing-python
* Installing conda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
* Installing Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
	+ Note, if you are on a remote server then you might not have permissions to create new conda environments where conda/mamba tries to store them. So you'll need to create a new place to store these enviroments and then when creating an environment suitable for snakemake use the `-p` parameter - e.g. `mamba create -p MY_ENV_PATH -c conda-forge -c bioconda -n snakemake snakemake`
* Installing tmux: https://github.com/tmux/tmux/wiki

Before starting any analyses, you should install the required R packages by running the [`install-packages.R`](install-packages.R) script, e.g. 

``` bash
Rscript install-packages.R
```

## To note

* The scripts will assume certain formats for various files. For example, it is assumed the cleaned version of your DNAm matrix will be an R object file (e.g. a ".RData" file) that can be loaded into R using the `load` function.
* The cell-type specific methods require that all covariates are numeric. We expect for the variable indicating AD case-control status, cases = 1 and controls = 0. For sex, we expect female = 1 and male = 2.

## Uploading the data

Please zip and compress your working directory and send it back. Command to zip/compress:

```bash 
tar -zcvf file.tar.gz /path/to/dir/
```

How to send the results can be found in the analysis plan you should have already received.