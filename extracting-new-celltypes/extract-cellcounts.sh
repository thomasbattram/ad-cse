#!/bin/bash

#SBATCH --job-name=extract_cellcounts
#SBATCH --account=smed001801
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=32GB

## Set working directory
wdir="" # FILL THIS IN
scriptdir="" # FILL THIS IN
cd ${wdir}

# url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167998/suppl/GSE167998_RAW.tar"
aries_dir="" # FILL THIS IN
tp="15up" # UPDATE ME
mkdir -p ${tp}
samp="${tp}/samplesheet.csv"
qc="${tp}/qc.objects.Robj"
cc="../data-extraction-and-qc/data/extended-blood-celltypes-epic-${tp}-idats.tsv"
idats="MA028" # UPDATE ME

Rscript ${scriptdir}/extract-cellcounts.R ${aries_dir} ${samp} ${qc} ${cc} ${idats} ${tp}