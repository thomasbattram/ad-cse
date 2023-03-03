# -------------------------------------------------------------------------------------
# Check AD definition associates with FLG
# ------------------------------------------------------------------------------------- 

## Null mutations in FLG are strong risk factors for AD, in this script we check whether our definition
## of AD associates strongly with FLG to see whether the definition is likely capturing AD cases

# srun --job-name "InteractiveJob" --account=smed001801 --partition=test --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=4GB --pty bash

## pkgs
library(aries) # get aries samples data
library(tidyverse) # tidy code and data
library(haven) # read in dta files
library(usefunc) # personal package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_samples <- args[2]
aries_dir <- args[3]
flg_file <- args[4]
outfile <- args[5]

## manual args
phen_file <- "data/ad-data.tsv"
meth_samples <- "data/samplenames.txt"
aries_dir <- ""
flg_file <- "data/children_FLG_variables.dta"
outfile <- "data/flg-associations.tsv"

## data
phen <- read_tsv(phen_file)
samples <- readLines(meth_samples)
flg_vars <- read_dta(flg_file)
# flg_vars <- flg_vars[, c("aln", "qlet", "FLG_comb")]

#
#
#

aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
samplesheet <- aries$samples %>%
	dplyr::filter(Sample_Name %in% samples)

all_phen <- phen %>%
	left_join(flg_vars) %>%
	left_join(samplesheet) %>%
	dplyr::select(Sample_Name, aln, alnqlet, qlet, ad, FLG_comb) %>%
	dplyr::filter(!is.na(ad), !is.na(FLG_comb))


#' Extract res from glm() function
#' 
#' @param glm_obj object obtained from running the glm() function
#' @return table containing summary stats
summ_glm_res <- function(glm_obj)
{
	summ_res <- summary(glm_obj)
	out <- tibble(Beta = summ_res$coef[2, 1], 
				  SE = summ_res$coef[2, 2], 
				  P = summ_res$coef[2, 4], 
				  OR = exp(Beta))
	return(out)
}


form <- as.formula(paste0("ad", " ~ ", "FLG_comb"))
glm_res <- glm(form, data = all_phen, family = "binomial")
summ_stats <- summ_glm_res(glm_res)
out <- summ_stats %>%
	mutate(N = nrow(all_phen), N_cases = sum(all_phen[["ad"]]), N_controls = N - N_cases)

write.table(out, file = outfile, col.names = T, row.names = F, sep = "\t", quote = F)