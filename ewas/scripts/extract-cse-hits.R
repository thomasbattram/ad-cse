# ----------------------------------------
# Extract hits from cell-specific effects methods
# ----------------------------------------

## Aim: Take EWAS results and extract any associations that replicate across methods

## Date: 2022-08-12

## pkgs
library(tidyverse) # tidy code, data, plots
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
celldmc_file <- args[1]
tca_file <- args[2]
outfile <- args[3]

## args testers
# celldmc_file <- "results/ewas/celldmc-res.RData"
# tca_file <- "results/ewas/tca-res.RData"
# outfile <- "results/celldmc-tca-hits.RData"

## data
celldmc_res <- new_load(celldmc_file)
tca_res <- new_load(tca_file)

# ----------------------------------------
# Setup
# ----------------------------------------

str(celldmc_res)
str(tca_res)

celltypes <- colnames(celldmc_res$beta)
#' Combine the summary statistics from each cell type
#' 
#' @param res results from the EWAS
#' @param celltype the cell type to combine results for
#' 
#' @return tibble of results with the summary statistics and "CpG" as columns 
comb_res <- function(res, celltype) 
{
	new_res <- map(res, as_tibble, rownames = "CpG")
	out <- map_dfc(new_res, celltype) %>%
		mutate(CpG = new_res$p$CpG) %>%
		dplyr::select(CpG, everything())
	return(out)
}

# lapply(celltypes, comb_res, res = tca_res)

#' Extract hits from EWAS res
#' 
#' @param res results from the EWAS
#' @param p_threshold P value threshold to extract associations
#' @param celltype the cell type to combine results for
#' 
#' @return tibble of results with the summary statistics
extract_hits <- function(res, p_threshold, celltype)
{
	tab_res <- comb_res(res, celltype)
	tab_res %>%
		dplyr::filter(p < p_threshold)	
}

#' Extract CpGs to check for replication from EWAS res
#' 
#' @param res results from the EWAS
#' @param cpgs CpGs of interest
#' @param celltype the cell type to combine results for
#' 
#' @return tibble of results with the summary statistics
extract_rep <- function(res, cpgs, celltype)
{
	tab_res <- comb_res(res, celltype)
	tab_res %>%
		dplyr::filter(CpG %in% cpgs)
}

# ----------------------------------------
# Extract the hits
# ----------------------------------------

## Look at replication in both directions, i.e.
## 1. Extract hits from CellDMC results
## 2. Check for replication in TCA results
## 3. Repeat steps reversing the methods


disc_method <- c("celldmc", "tca")
rep_method <- c("tca", "celldmc")
res_to_extract <- expand.grid(disc = disc_method, 
							  rep = rep_method, 
							  celltype = celltypes, 
							  stringsAsFactors = FALSE)

res_to_extract <- res_to_extract %>%
	dplyr::filter(disc != rep)

n_celltypes <- length(celltypes)
p_threshold <- 1e-5

all_res <- lapply(1:nrow(res_to_extract), function(x) {
	temp <- res_to_extract[x, ]
	## extract hits from discovery res
	disc_res <- get(paste0(temp$disc, "_res"))
	disc_hits <- extract_hits(disc_res, p_threshold, temp$celltype)
	if (nrow(disc_hits) == 0) return(NULL)
	## Extract hits from replication results and check for replication
	rep_res <- get(paste0(temp$rep, "_res"))
	rep_hits <- extract_rep(rep_res, cpgs = disc_hits$CpG, temp$celltype)
	comb_hits <- disc_hits %>%
		left_join(rep_hits, by = c("CpG"), suffix = c("_disc", "_rep")) %>%
		mutate(p_rep_fdr = p.adjust(p_rep, n = n_celltypes * nrow(disc_hits), method = "fdr")) %>%
		mutate(replicated = ifelse(p_rep_fdr < 0.05 & sign(beta_disc) == sign(beta_rep), TRUE, FALSE))
	return(comb_hits)
})

names(all_res) <- paste0(res_to_extract$disc, "_", res_to_extract$rep, "_", res_to_extract$celltype)

## Quick check of data
all_res[!sapply(all_res, is.null)]

out <- list(res = all_res[!sapply(all_res, is.null)], 
			note = "names of results are in the format: discovery-method_replication-method_celltype")

save(out, file = outfile)


