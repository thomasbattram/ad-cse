# ------------------------------------------------
# Extract ARIES samples
# ------------------------------------------------

## Aim: Extract ARIES samples from 15up data

## This needs to be run on either BC4 or BP1

## pkgs
library(aries)

## args
aries_dir <- "" ## FILL THIS IN 
outfile_15up <- "data/aries-samples-15up.tsv"
outfile_f7 <- "data/aries-samples-f7.tsv"

## Get data
aries <- aries.select(aries_dir, time.point="15up", featureset = "epic")
samples <- aries$samples[, c("aln", "qlet", "Sample_Name")]

write.table(samples, outfile_15up, col.names = T, row.names = F, quote = F, sep = "\t")

aries.time.points(aries_dir)
aries <- aries.select(aries_dir, time.point="F7")
samples <- aries$samples[, c("aln", "qlet", "Sample_Name")]

write.table(samples, outfile_f7, col.names = T, row.names = F, quote = F, sep = "\t")