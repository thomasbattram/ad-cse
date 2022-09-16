# ---------------------------------------------------------------
# Estimate extended blood cell counts 
# ---------------------------------------------------------------

## Aim: Estimate whole blood cell counts in 15up individuals using the new reference with 12 cell types

## NOTE: Code comes from: https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r

# srun --job-name "InteractiveJob" --partition=veryshort --nodes=1 --ntasks-per-node=4 --cpus-per-task=4 --time=6:00:00 --mem=50GB --pty bash

## pkgs
library(tidyverse)
library(aries)
library(meffil)
library(minfi)
library(GEOquery)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## args
args <- commandArgs(trailingOnly = TRUE)
acc <- args[1]
url <- args[2]
aries_dir <- args[3]
geo_outfile <- args[4]
outfile <- args[5]
minfi_functions <- args[6]

## Testing args
# acc <- "GSE167998"
# url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167998/suppl/GSE167998_RAW.tar"
# aries_dir <- ""
# geo_outfile <- "GSE167998/samples.csv"
# outfile <- "../data-extraction-and-qc/extended-blood-celltypes-epic-15up.tsv"
# minfi_functions <- "minfiEPIC.r"

source(minfi_functions)

dir.create(path <- acc)
if(length(list.files(path, "*.idat$")) == 0) {
    filename <-  file.path(path, basename(url))
    download.file(url, filename)
    cat(date(), "Extracting files from GEO archive.\n")
    system(paste("cd", path, ";", "tar xvf", basename(filename)))
    unlink(filename)
    cat(date(), "Unzipping IDAT files.\n")
    system(paste("cd", path, ";", "gunzip *.idat.gz"))

    geo <- getGEO(acc, GSEMatrix=F)
    geo <- lapply(geo@gsms, function(gsm) unlist(gsm@header))
    geo <- do.call(rbind, geo)
    geo <- as.data.frame(geo, stringAsFactors=F)
    write.csv(geo, file=geo_outfile)
}

samplesheet <- meffil::meffil.create.samplesheet(path)
geo_samples <- read.csv(geo_outfile)

# geo_samples[, c("X", "data_row_count")]

celltypes <- c("Neu", "Eos", "Bas", "Mono", 
               "Bnv", "Bmem", "CD4nv", "CD4mem", 
               "Treg", "CD8nv", "CD8mem", "NK", "MIX")

samples <- data.frame(Sample_Name = geo_samples$X, CellType = NA)
for (i in 1:nrow(geo_samples)) {
    ct <- sapply(celltypes, function(x) {
        any(grepl(x, geo_samples[i, c("data_row_count", "description", "data_processing3")]))
    }) 
    samples[i, "CellType"] <- names(ct)[ct]
}
table(samples$CellType)

# qc.objects <- meffil::meffil.qc(samplesheet, cell.type.reference=NA, verbose=T)


## Extract data
RGset <- read.metharray.exp(targets = samplesheet, force=TRUE)

reference <- preprocessNoob(RGset)
reference <- mapToGenome(reference)

extracted.data <- extractFromRGSetEPIC(RGset)
samplesheet <- colData(RGset)
samplesheet <- merge(samplesheet, samples)
sex <- getSex(reference, cutoff = -3)$predictedSex
sex <- sign(sex=="F")+1
reference <- normalizeFunnormEPIC(
    object=reference, extractedData=extracted.data,
    sex=sex, nPCs=10, verbose=F)
M <- getMeth(reference)
U <- getUnmeth(reference)

## Add cell type reference to meffil
cell.types <- celltypes[celltypes != "MIX"]
stopifnot(all(cell.types %in% samplesheet$CellType))
selected <- samplesheet$CellType %in% cell.types
meffil::meffil.add.cell.type.reference(
    "blood extended idoloptimized epic",
    M[,selected], U[,selected],
    cell.types=samplesheet$CellType[selected],
    chip="epic",
    featureset="epic",
    # specific.sites=IDOLOptimizedCpGs,
    # description="Derived from FlowSorted.BloodExtended.EPIC",
    verbose=T)

meffil.list.cell.type.references()

## remove clutter
rm(list = c("RGset", "reference", "extracted.data", "samplesheet", "sex", "M", "U", "geo_samples", "selected"))

## Estimate cell counts in 15up individuals
aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
# test_samp <- head(aries$samples)
beta <- aries.methylation(aries)
meth <- beta[, aries$samples$Sample_Name]
rm(beta)

out_celltypes <- meffil.estimate.cell.counts.from.betas(meth, "blood extended idoloptimized epic")

out_celltypes <- as.data.frame(out_celltypes)
out_celltypes$Sample_Name <- rownames(out_celltypes)

## Write it out
write.table(out_celltypes, file = outfile, sep = "\t", row.names = F, col.names = T, quote = F)
