# -------------------------------------------------------
# Filter ARIES betas
# -------------------------------------------------------
# Version = v4

# -------------------------------------------------------
# Setup
# -------------------------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(meffil) # contains 450k features
library(usefunc) # own package of useful functions

# -------------------------------------------------------
# Load ARIES data and remove bad samples
# -------------------------------------------------------

## load aries data
samplesheet <- new_load(snakemake@input[["samplesheet_file"]])
samplesheet <- samplesheet %>%
	dplyr::filter(time_point == snakemake@params[["timepoint"]])

# samples to remove
sample_rm <- which(samplesheet$duplicate.rm == "Remove" | samplesheet$genotypeQCkids == "ETHNICITY" | 
				   samplesheet$genotypeQCkids == "HZT;ETHNICITY" | samplesheet$genotypeQCmums == "/strat")
samplesheet <- samplesheet[-sample_rm, ]

# methylation data 
beta <- new_load(snakemake@input[["aries_meth_dat"]])
meth <- beta[, samplesheet$Sample_Name]
rm(beta)

# detection p values
detp <- new_load(snakemake@input[["aries_detection_p"]])
pvals <- detp[, samplesheet$Sample_Name]
rm(detp)

print("finished reading in stuff")

# -------------------------------------------------------
# Remove bad probes
# -------------------------------------------------------

## load annotation data
annotation <- meffil.get.features("450k")

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

# -------------------------------------------------------
# Filter ARIES betas
# -------------------------------------------------------

q <- rowQuantiles(meth, probs = c(0.25, 0.75), na.rm = T)
iqr <- q[, 2] - q[, 1]
too.hi <- which(meth > q[,2] + 3 * iqr, arr.ind=T)
too.lo <- which(meth < q[,1] - 3 * iqr, arr.ind=T)
if (nrow(too.hi) > 0) meth[too.hi] <- NA
if (nrow(too.lo) > 0) meth[too.lo] <- NA

dim(meth)
num_na <- apply(meth, 2, function(x){sum(is.na(x))})
rem_samp <- which(num_na > (0.05 * nrow(meth)))
meth <- meth[, -rem_samp]
dim(meth)

print(paste0("Number of samples removed = ", length(rem_samp)))

save(meth, file = snakemake@output[[1]])

