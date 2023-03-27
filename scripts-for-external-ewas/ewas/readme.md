# EWAS

Aim: Run conventional EWAS, cell-specific EWAS, and variance EWAS of AD

## Before running any scripts

All phenotype and DNA methylation data should have been cleaned.

please make the required folders in your working directory

In linux:

``` bash
# Change WORK_DIR if you need to 
WORK_DIR="."
cd ${WORK_DIR}
mkdir -p data/svs
mkdir -p results/ewas/tca-temp
```

## Snakemake workflow

If you're using snakemake also run `mkdir job-reports`

The scripts can be executed manually or via snakemake. Here is the snakemake workflow:

0. Start a tmux session - `tmux new -s ad-cse`
1. Activate conda env - `conda activate snakemake`
2. Edit the [`Snakefile`](Snakefile)
3. Do a dry run with `snakemake -nrp`
4. If you're running locally, then run:

```bash
snakemake -rp --cores 1
```

If you're running on a high performance computing (HPC) environment then edit `cluster.json`, and edit the command the below to match your HPC specifications. The example below is for a `slurm` environment.

``` bash
snakemake -rp \
-j 1 \
--cluster-config cluster.json \
--cluster "sbatch \
        --job-name={cluster.name} \
        --partition={cluster.partition} \
        --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} \
        --cpus-per-task={cluster.ncpu} \
        --time={cluster.time} \
        --mem={cluster.mem} \
        --output={cluster.output} \
        --error={cluster.error}"
```

5. Deactivate the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t ad-cse`

## Script order

* [gen-svs.R](scripts/gen-svs.R) : Generate SVs 
* [conventional-ewas.R](scripts/conventional-ewas.R) : Run the EWAS using ewaff
* [celldmc-ewas.R](scripts/celldmc-ewas.R) : Run the EWAS using celldmc
* [tca-ewas.R](scripts/tca-ewas.R) : Run the EWAS using TCA (split from celldmc script because TCA requires splitting the DNAm data into chunks)
* [combine-tca-res.R](scripts/combine-tca-res.R): Combine TCA results into a single file
* [var-ewas.R](scripts/var-ewas.R) : Run variance EWAS

## Running without Snakemake

The scripts can be run without Snakemake. This can be done manually using R (note that some take many hours to run), or by running them using the command `RScript` on a linux system (which could be an HPC environment). There is an example of how this is done for the conventional EWAS in [`conv-ewas-example.sh`](scripts/conv-ewas-example.sh) with the equivalent part of the [`Snakefile`](Snakefile) for comparison. 

The scripts should be run in the order shown above under the heading `Script order`

If running without Snakemake it is important to remember we want to run the following EWAS:

1. Conventional EWAS adjusted for cell counts and unadjusted for cell counts. Outputs = results/ewas/ewaff-res-cc.tsv AND results/ewas/ewaff-res-no-cc.tsv.
2. CellDMC EWAS. Output = results/ewas/celldmc-res.RData
3. TCA EWAS. Output = results/ewas/tca-res.RData
4. Variance EWAS adjusted for cell counts, unadjusted for cell counts, and with no covariates. Outputs = results/ewas/varewas-res-cc.tsv results/ewas/varewas-res-no-cc.tsv results/ewas/varewas-res-no-covs.tsv



