# EWAS

Aim: Run conventional and cell-specific EWAS of AD

Requires having run the scripts in the [data-extraction-and-qc](../[data-extraction-and-qc]) folder.


* [gen-svs.R](scripts/gen-svs.R) : Generate SVs 
* [conventional-ewas.R](scripts/conventional-ewas.R) : Run the EWAS using ewaff
* [cse-ewas.R](scripts/cse-ewas.R) : Run the EWAS using cell-dmc and TCA
* [extract-cse-hits.R](scripts/extract-cse-hits.R) : Extract all sites that are associated with AD in a cell-type with one method and replicate in the other method
* [omicwas-ewas.R](scripts/omicwas-ewas.R) : Run omicWAS with all the replicated sites
* [summarise-results.R](scripts/summarise-results.R) : Analyse the results (Manhattans, etc.) + output nice tables 
* [ewas-report.Rmd](report/ewas-report.Rmd) : Put results into a report


## NOTE: NEED TO ANNOTATE SOME SCRIPTS PROPERLY!!