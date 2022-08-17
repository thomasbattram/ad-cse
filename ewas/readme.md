# EWAS

Aim: Run conventional and cell-specific EWAS of AD

Requires having run the scripts in the [data-extraction-and-qc](../[data-extraction-and-qc]) folder.

The scripts can be executed manually or via snakemake. Here is the snakemake workflow:

0. Start a tmux session - `tmux new -s ad-cse`
1. Activate conda env - `conda activate /user/home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:

``` bash
snakemake -p \
-j 1 \
--cluster-config bc4-cluster.json \
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

The scripts run by the snakemake workflow are below:

* [gen-svs.R](scripts/gen-svs.R) : Generate SVs 
* [conventional-ewas.R](scripts/conventional-ewas.R) : Run the EWAS using ewaff
* [cse-ewas.R](scripts/cse-ewas.R) : Run the EWAS using cell-dmc and TCA
* [extract-cse-hits.R](scripts/extract-cse-hits.R) : Extract all sites that are associated with AD in a cell-type with one method and replicate in the other method
* [omicwas-ewas.R](scripts/omicwas-ewas.R) : Run omicWAS with all the replicated sites
* [summarise-results.R](scripts/summarise-results.R) : Analyse the results (Manhattans, etc.) + output nice tables 
* [ewas-report.Rmd](report/ewas-report.Rmd) : Put results into a report