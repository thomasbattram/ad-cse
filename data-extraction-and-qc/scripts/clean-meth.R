
# -------------------------------------------------------
# Filter ARIES betas
# -------------------------------------------------------
# Version = v4

# srun --job-name "InteractiveJob" --partition=veryshort --nodes=1 --ntasks-per-node=4 --cpus-per-task=4 --time=6:00:00 --mem=50GB --pty bash

# -------------------------------------------------------
# Setup
# -------------------------------------------------------

## pkgs
library(tidyverse) # tidy code and data
# library(meffil) # contains DNAm data annotations - doesn't work on bc4 because Cairo isn't installed... - 
library(aries) # extract aries data
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
aries_dir <- args[2]
zhou_file <- args[3]
outfile_meth <- args[4]
outfile_samples <- args[5]

## args testers
# phen_file <- "data/ad-data.tsv"
# aries_dir <- "/user/work/ms13525/aries"
# zhou_file <- "data/epic-zhou-list.txt"
# outfile_meth <- "data/clean-meth.RData" 
# outfile_samples <- "data/samplenames.txt"

# -------------------------------------------------------
# Load ARIES data and remove bad samples
# -------------------------------------------------------

## load phenotype file
pheno <- read_tsv(phen_file)

## load aries data
aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")

samplesheet <- aries$samples %>%
	dplyr::filter(time_point == "15up") %>%
	left_join(pheno) %>%
	dplyr::filter(!is.na(ad)) %>%
	dplyr::select(-ad)

rm(pheno)

# samples to remove - need to find out how to spot sex mismatches...
sample_rm <- which(samplesheet$duplicate.rm | samplesheet$genotype.mismatch)
samplesheet <- samplesheet[-sample_rm, ]

# methylation data 
beta <- aries.methylation(aries)
meth <- beta[, samplesheet$Sample_Name]
rm(beta)

# detection p values
detp <- aries.detectionp(aries)
pvals <- detp[, samplesheet$Sample_Name]
rm(detp)

message("finished reading in aries DNAm data")

# -------------------------------------------------------
# Remove bad probes
# -------------------------------------------------------

## load annotation data
annotation <- meffil::meffil.get.features("epic")

## Filter meth data (remove sex chromosomes and SNPs and probes with high detection P-values)
pvalue_over_0.05 <- pvals > 0.05
count_over_0.05 <- rowSums(sign(pvalue_over_0.05))
Probes_to_exclude_Pvalue <- rownames(pvals)[which(count_over_0.05 > ncol(pvals) * 0.05)]
XY <- as.character(annotation$name[which(annotation$chromosome %in% c("chrX", "chrY"))])
SNPs.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])
annotation<- annotation[-which(annotation$name %in% c(XY, SNPs.and.controls, Probes_to_exclude_Pvalue)), ]
print(length(annotation))
print(dim(meth))
meth <- base::subset(meth, row.names(meth) %in% annotation$name)
paste("There are now ", nrow(meth), " probes")
paste(length(XY), "were removed because they were XY")
paste(length(SNPs.and.controls), "were removed because they were SNPs/controls")
paste(length(Probes_to_exclude_Pvalue), "were removed because they had a high detection P-value")
rm(XY, SNPs.and.controls, pvals, count_over_0.05, pvalue_over_0.05, Probes_to_exclude_Pvalue)

filtered_vars <- c("detection_p_values", "on_XY", "SNPs/controls")

# COULD ALSO ADD ZHOU LIST HERE! 
epic_zhou <- usefunc::get_zhou_recs(outpath = zhou_file, array = "epic")
meth <- meth[!rownames(meth) %in% epic_zhou, ]

# -------------------------------------------------------
# Filter ARIES samples
# -------------------------------------------------------

dim(meth)
num_na <- apply(meth, 2, function(x){sum(is.na(x))})
rem_samp <- which(num_na > (0.05 * nrow(meth)))
meth <- meth[, -rem_samp]
dim(meth)

print(paste0("Number of samples removed = ", length(rem_samp)))

save(meth, file = outfile_meth)
writeLines(colnames(meth), con = outfile_samples)
