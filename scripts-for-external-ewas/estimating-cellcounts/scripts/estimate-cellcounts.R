# ---------------------------------------------------------------
# Estimate extended blood cell counts 
# ---------------------------------------------------------------

## Aim: Estimate counts from 12 cell types using whole blood DNAm data

## NOTE: Code comes from: https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r

## pkgs
library(meffil)

## args
args <- commandArgs(trailingOnly = TRUE)
samp_outfile <- args[1] # samplesheet output filename
qc_outfile <- args[2] # meffil QC objects output filename
cc_outfile <- args[3] # cell count output filename
path_to_idat_files <- args[4] # path to directory containing IDAT files

## manual args
# samp_outfile <- "data/samplesheet.csv" # samplesheet output filename
# qc_outfile <- "data/qc.objects.Robj" # meffil QC objects output filename
# cc_outfile <- "data/extended-blood-celltypes-idats.tsv" # cell count output filename
# path_to_idat_files <- "" # path to directory containing IDAT files

samplesheet <- meffil.create.samplesheet(path_to_idat_files, recursive=TRUE)
write.csv(samplesheet, file = samp_outfile, quote = TRUE, row.names = FALSE)

## Checking the correct reference is stored in meffil
if (!"blood gse167998" %in% meffil.list.cell.type.references()) {
	stop("The required reference is not in the version of meffil you are using, please update it using:\n", 
		 "devtools::install_github('perishky/meffil')")
}

qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse167998", verbose=TRUE)
save(qc.objects, file= qc_outfile)

## Check to see whether any of the QC objects were errors
# which(sapply(qc.objects, class) == 'try-error')

# load(qc_outfile)

## Output the cell counts
cc <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc <- data.frame(IID = row.names(cc), cc)
write.table(cc, file = cc_outfile, sep = "\t", row.names = F, col.names = T, quote = F)