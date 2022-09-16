#!/bin/bash

#SBATCH --job-name=extract_cellcounts
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --mem=50GB

## Set working directory
wdir=""
cd ${wdir}

acc="GSE167998"
url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167998/suppl/GSE167998_RAW.tar"
aries_dir="" # FILL THIS IN
geo_outfile="GSE167998/samples.csv"
outfile="../data-extraction-and-qc/extended-blood-celltypes-epic-15up.tsv"
minfi_functions="minfiEPIC.r"


Rscript extract-cellcounts.R ${acc} ${url} ${aries_dir} ${geo_outfile} ${outfile} ${minfi_functions}