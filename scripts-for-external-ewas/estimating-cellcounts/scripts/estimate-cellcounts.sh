#!/bin/bash

#SBATCH --job-name=extract_cellcounts
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

mkdir -p "data"
samp="data/samplesheet.csv"
qc="data/qc.objects.Robj"
cc="data/extended-blood-celltypes-epic-idats.tsv"
idats="" # FILL THIS IN

Rscript ${scriptdir}/extract-cellcounts.R ${samp} ${qc} ${cc} ${idats}