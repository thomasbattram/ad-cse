# ---------------------------------------------------------------
# Estimate extended blood cell counts 
# ---------------------------------------------------------------

## Aim: Estimate whole blood cell counts in 15up individuals using the new reference with 12 cell types

## NOTE: Code comes from: https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r

## Cell counts are now in meffil


## pkgs
library(tidyverse)
library(aries)
library(meffil)

## args
args <- commandArgs(trailingOnly = TRUE)
aries_dir <- args[1]
samp_outfile <- args[2]
qc_outfile <- args[3]
cc_outfile <- args[4]
path_to_idat_files <- args[5]
# rsync -rv /projects/MRC-IEU/research/data/alspac/epigenetic/methylation/aries/raw/MA028/ MA028
aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
aries_samplesheet <- aries$samples
samplesheet <- meffil.create.samplesheet(path_to_idat_files, recursive=TRUE)
all(aries_samplesheet$Sample_Name %in% samplesheet$Sample_Name) # TRUE
write.csv(samplesheet, file = samp_outfile, quote = TRUE, row.names = FALSE)

## Checking the correct reference is stored in meffil
"blood gse167998" %in% meffil.list.cell.type.references()

# which(sapply(qc.objects, class) == 'try-error')[1]
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
save(qc.objects, file= qc_outfile)

# load(qc_outfile)

cc <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc <- data.frame(Sample_Name = row.names(cc), cc)
write.table(cc, file = cc_outfile, sep = "\t", row.names = F, col.names = T, quote = F)

# ## Estimate cell counts in 15up individuals
# aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
# # test_samp <- head(aries$samples)
# # beta <- aries.methylation(aries)
# # meth <- beta[, test_samp$Sample_Name]
# beta <- aries.methylation(aries)
# meth <- beta[, aries$samples$Sample_Name]
# rm(beta)

# out_celltypes <- meffil.estimate.cell.counts.from.betas(meth, "blood extended idoloptimized epic")

# out_celltypes <- as.data.frame(out_celltypes)
# out_celltypes$Sample_Name <- rownames(out_celltypes)

# ## Write it out
# write.table(out_celltypes, file = outfile, sep = "\t", row.names = F, col.names = T, quote = F)
