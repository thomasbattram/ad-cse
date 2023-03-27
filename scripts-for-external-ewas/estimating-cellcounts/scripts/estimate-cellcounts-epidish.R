# ---------------------------------------------------------------
# Estimate extended blood cell counts 
# ---------------------------------------------------------------

## Aim: Estimate counts from 12 cell types using whole blood DNAm data

## NOTE: Code comes from: https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r

## pkgs
library(tidyverse)
library(matrixStats) 
library(EpiDISH)
library(usefunc)

## args
args <- commandArgs(trailingOnly = TRUE)
dnam_mat_file <- args[1] # DNAm matrix file
cc_outfile <- args[2] # cell count output filename

## manual args
# dnam_mat_file <- "data/clean-meth.RData" # DNAm matrix file
# cc_outfile <- "data/extended-blood-celltypes-epidish.tsv" # cell count output filename

## data
data(cent12CT.m)
meth <- new_load(dnam_mat_file)

## Estimate cell counts
meth <- impute_matrix(meth)
out.l <- epidish(beta.m = meth, ref.m = cent12CT.m, method = "RPC")
epi_rpc <- as_tibble(out.l$estF, rownames = "IID") %>%
	rename(CD4nv = CD4Tnv, Bas = Baso, CD4mem = CD4Tmem, CD8mem = CD8Tmem, 
		   CD8nv = CD8Tnv)

write.table(epi_rpc, cc_outfile, row.names = F, col.names = T, quote = F, sep = "\t")