## From https://github.com/perishky/meffil/blob/master/data-raw/idoloptimized-references.r

#### TESTING WHAT IT LOOKS LIKE FOR OTHER REFERENCES

library(meffil)
require("minfi")
library(GEOquery)
source("minfiEPIC.r")
require("FlowSorted.Blood.EPIC")

## module add languages/r/4.1.0
require("ExperimentHub")
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
ds <- hub[["EH1136"]]

## WHAT DOES DS LOOK LIKE?!?! - rgset - need to use minfi to run it all
str(ds)

reference <- preprocessNoob(ds)
reference <- mapToGenome(reference)

#### TEST OVER
