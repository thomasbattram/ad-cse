# omicWAS replication

There is one script in this folder, and that runs a cell-type specific EWAS using omicWAS on any associations found within the meta-analysis of results from cell-type specific EWAS using CellDMC and TCA. 

## Before running any scripts

All phenotype and DNA methylation data should have been cleaned (use the same data as from the previous round of analyses). We assume your working directory setup is the same as before.

You should have been sent the file file `meta-analysis-hits.txt`. Please move this into the `data` folder in your working directory.

## Workflow

Run the R code manually or you can run the following bash command

```bash
Rscript scripts/omicwas-ewas.R "data/cleaned-pheno-data.tsv" \
							   "data/clean-meth.RData" \
							   "data/svs/ad-svs.tsv" \
							   "data/covars-no-cc.txt" \
							   "data/meta-analysis-hits.txt" \
							   "results/ewas/omicwas-res.RData"
```

