# ------------------------------------------------
# Extract ARIES samples
# ------------------------------------------------

## Aim: Extract ARIES samples from 15up data

## This needs to be run on either BC4 or BP1

## pkgs
library(aries)

## args
aries_dir <- "" ## FILL THIS IN 
outfile <- "data/aries-samples-15up.tsv"

## Get data
aries <- aries.select(aries_dir, time.point="15up", featureset = "epic")
samples <- aries$samples[, c("aln", "qlet", "Sample_Name")]

write.table(samples, outfile, col.names = T, row.names = F, quote = F, sep = "\t")