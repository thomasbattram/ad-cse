# ------------------------------------------------------
# Testing celldmc res
# ------------------------------------------------------

## Date: 2022-08-21

## pkgs
library(tidyverse)
library(usefunc)
library(meta) # for meta-analysis


## args
args <- commandArgs(trailingOnly = TRUE)


## args tester
# celldmc_file <- "results/ewas/celldmc-res.RData"
# ewaff_cc_file <- "results/ewas/ewaff-res-cc.tsv"
# ewaff_no_cc_file <- "results/ewas/ewaff-res-no-cc.tsv"


## data
celldmc <- new_load(celldmc_file)
ewaff_cc <- read_tsv(ewaff_cc_file)
ewaff_no_cc <- read_tsv(ewaff_no_cc_file)

ewaff_cc <- ewaff_cc %>%
	mutate(celltype = "ewaff_cc") %>%
	dplyr::select(CpG = probeID, celltype, beta = BETA, se = SE, p = P)
ewaff_no_cc <- ewaff_no_cc %>%
	mutate(celltype = "ewaff_no_cc") %>%
	dplyr::select(CpG = probeID, celltype, beta = BETA, se = SE, p = P)

## 
comb_res <- function(res, celltype) {
	new_res <- map(res, as_tibble, rownames = "CpG")
	out <- map_dfc(new_res, celltype) %>%
		mutate(CpG = new_res$p$CpG) %>%
		dplyr::select(CpG, everything())
	return(out)
}

celltypes <- names(celldmc$beta)
celldmc_res <- lapply(celltypes, comb_res, res = celldmc)
names(celldmc_res) <- celltypes
celldmc_res <- bind_rows(celldmc_res, .id = "celltype")
cpgs <- unique(celldmc_res$CpG)
set.seed(2)
samp_cpgs <- sample(cpgs, 1000)

start_time <- proc.time()
res <- lapply(samp_cpgs, function(cpg) {
	print(cpg)
	celldmc_res_filt <- celldmc_res[celldmc_res$CpG %in% cpg, ]
	m_gen <- metagen(TE = beta,
					 seTE = se,
					 studlab = celltype,
					 data = celldmc_res_filt,
					 sm = "SMD",
					 fixed = FALSE,
					 random = TRUE,
					 method.tau = "REML",
					 hakn = TRUE,
					 title = "")
	meta_out <- tibble(celltype = c("meta-fixed", "meta-random"),
					   CpG = c(cpg, cpg),
					   beta = c(m_gen$TE.fixed, m_gen$TE.random), 
					   se = c(m_gen$seTE.fixed, m_gen$seTE.random),
					   p = c(m_gen$pval.fixed, m_gen$pval.random))
	ewaff_out <- ewaff_cc %>%
		dplyr::filter(CpG == cpg) %>%
		bind_rows(ewaff_no_cc[ewaff_no_cc$CpG == cpg,])
	out <- bind_rows(celldmc_res_filt, meta_out, ewaff_out)
	return(out)
})
time_taken <- proc.time() - start_time
time_taken # elappsed = 253.992 (~4mins)


all_res <- bind_rows(res) %>%
	dplyr::filter(celltype %in% c("meta-fixed", "meta-random", "ewaff_cc", "ewaff_no_cc"))

methods <- unique(all_res$celltype)
comp_res <- map_dfc(methods, function(method) {
	all_res %>%
		dplyr::filter(celltype == method) %>%
		pull(beta)
})
colnames(comp_res) <- methods

cor(comp_res)