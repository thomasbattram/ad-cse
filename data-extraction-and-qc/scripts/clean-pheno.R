# -------------------------------------------------
# Clean phenotype data for EWAS
# -------------------------------------------------

## Aim: take phenotype data and clean it ready for EWAS

## pkgs
library(tidyverse) # tidy data, code, plots
library(aries) # aries samplesheet
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
phen_infile <- args[1]
meth_samples <- args[2]
cellcount_file <- args[3] # keep this as "" if using "aries" package
aries_dir <- args[4]
phen_outfile <- args[5]

## args testers
# phen_infile <- "data/ad-data.tsv"
# meth_samples <- "data/samplenames.txt"
# cellcount_file <- "data/extended-blood-celltypes-epic-15up.tsv"
# aries_dir <- "/user/work/ms13525/aries"
# phen_outfile <- "data/ad-data-cleaned.tsv"

## data
samples <- readLines(meth_samples)
phen <- read_tsv(phen_infile)

# -------------------------------------------------------
# Load ARIES samplesheet and extract only samples in DNAm data
# -------------------------------------------------------

aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
samplesheet <- aries$samples %>%
	dplyr::filter(Sample_Name %in% samples)

if (cellcount_file == "") {
	cell_counts <- aries$cell.counts[["blood-gse35069-complete"]] %>%
		as.data.frame() %>%
		rownames_to_column(var = "Sample_Name") %>%
		as_tibble()
} else {
	cell_counts <- read_tsv(cellcount_file)
}

pheno <- phen %>%
	left_join(samplesheet) %>%
	left_join(cell_counts) %>%
	dplyr::filter(!is.na(Sample_Name)) %>%
	dplyr::select(Sample_Name, aln, alnqlet, qlet, ad, age, sex, all_of(colnames(cell_counts)))

write.table(pheno, file = phen_outfile, quote = F, col.names = T, row.names = F, sep = "\t")