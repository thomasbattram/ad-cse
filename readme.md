# Cell-specific associations between DNAm and AD

## For collaborators

For external collaborators, please clone this repository and then use [scripts-for-external-ewas](scripts-for-external-ewas) for all scripts and explanations of how these scripts work.

## Brief intro

Briefly, this project will include both conventional and cell-specific EWAS of AD in ARIES participants and potentially other cohorts. 

Analyses will be split into different folders, each with their own scripts and workflows. These will be detailed within readmes within those folders. 

## Analyses

* [extracting-new-celltypes](extracting-new-celltypes)
	+ Extract cell types from the Salas 2022 paper

* [data-extraction-and-qc](data-extraction-and-qc)
	+ Define AD in ARIES participants and QC pheno + DNAm data

* [ewas](ewas)
	+ Run conventional and cell-specific EWAS of AD

EXTRA STUFF!


## Notes

Initially the plan was to adjust for batch by generating SVs and adding them as covariates in the model. HOWEVER, this was problematic because the SVs were capturing cell counts also. So instead, batch effects were adjusted for directly by adjusting for "plate". Slide was considered when thinking of model adjustment, but DNAm PCs associate with slide and plate with similar strength, suggesting the variance captured by plate will include that captured by slide. 

Update AD-CSE report to reflect latest things of interest:
-	Add in cell count prediction differences in cases and controls (violin/box plot)
-	Highlight CpG sites that replicated in omicWAS to Manhattans
-	Update report to reflect new methods and results
	+	Also update plots so they don’t look shite
-	Make sure the adjustment for multiple testing includes 12 cell counts not 7
-	Maybe change identification of sites to P value threshold that isn’t 1x10-7!!!
