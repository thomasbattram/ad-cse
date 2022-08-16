# -------------------------------------------------------
# Script to generate surrogate variables for ARIES
# -------------------------------------------------------

# srun --job-name "InteractiveJob" --partition=veryshort --nodes=1 --ntasks-per-node=4 --cpus-per-task=4 --time=6:00:00 --mem=50GB --pty bash

## pkgs
library(tidyverse) # tidy data and code
library(sva) # calculating surrogate variables
library(SmartSVA) # calculating SVs
library(matrixStats) # imputing DNAm data
library(usefunc) # personal package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
out_file <- args[3]
removed_out_file <- args[4]

# phen_file <- "data/ad-data-cleaned.tsv"
# meth_file <- "data/clean-meth.RData"
# out_file <- "data/svs/ad-svs.tsv"
# removed_out_file <- "data/svs/ad-removed-svs.RData"

## read in data
phen_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)

# -------------------------------------------------------
# functions
# -------------------------------------------------------

#' Impute missing values in DNAm matrix
#' 
#' @param x DNAm matrix
#' @param FUN function to apply to "x" to get values to impute missing values
#' 
#' @return imputed DNAm matrix
impute_matrix <- function(x, FUN = function(x) rowMedians(x, na.rm = T)) {
    idx <- which(is.na(x), arr.ind = T)
    if (length(idx) > 0) {
        v <- FUN(x)
        v[which(is.na(v))] <- FUN(matrix(v, nrow = 1))
        x[idx] <- v[idx[, "row"]]
    }
    return(x)
}

# function to add quotes for weird trait names
addq <- function(x) paste0("`", x, "`")

#' Generate surrogate variables
#' 
#' @param trait trait of interest
#' @param phen_data phenotype data (data.frame or tibble) containing trait of interest and covariates
#' @param meth_data DNAm matrix
#' @param covariates character vector of covariates
#' @param nsv number of surrogate variables to use
#' @param IID name of the identifier used for DNAm samples. Default = "Sample_Name"
#' 
#' @return data.frame of surrogate variable values for each individual
generate_svs <- function(trait, phen_data, meth_data, covariates = "", nsv, 
						 IID = "Sample_Name") {
	print("Starting SV generation")
	if (any(grepl("^sv", paste(covariates, collapse = "|")))) {
		covs <- covariates[-grep("^sv", covariates)]
	} else {
		covs <- covariates
	}
	phen <- phen_data %>%
		dplyr::select(one_of(IID), one_of(trait, covs)) %>%
		.[complete.cases(.), ]
	
	mdat <- meth_data[, colnames(meth_data) %in% phen[[IID]]]
	phen <- phen %>%
		dplyr::filter(!!as.symbol(IID) %in% colnames(mdat))
	
	# models 
	trait_mod <- paste0("~ ", addq(trait))
	cov_mod <- paste(covs, collapse = " + ")
	if (covs != "") {
		full_mod <- paste(trait_mod, cov_mod, sep = " + ")
		fom <- as.formula(full_mod)
		# null model
		fom0 <- as.formula(paste0("~ ", cov_mod))
		mod0 <- model.matrix(fom0, data = phen)
	} else {
		fom <- as.formula(trait_mod)
		mod0 <- NULL
	}

	# full model - with variables of interest 
	mod <- model.matrix(fom, data = phen)

	# Estimate the surrogate variables
	tryCatch({
		svobj <- smartsva.cpp(mdat, mod, mod0, n.sv = nsv, VERBOSE = T)
		svs <- as.data.frame(svobj$sv, stringsAsFactors = F)
		svs[[IID]] <- phen[[IID]]
		# head(svs)
		colnames(svs)[1:nsv] <- paste0("sv", 1:nsv)
		return(svs)
	}, error = function(e) {err_msg(e, r_msg = TRUE, user_msg = paste("SV fail"))})
}

# -------------------------------------------------------
# sort data for generating SVs
# -------------------------------------------------------

# methylation data
mdata <- impute_matrix(meth)
rm(meth)

# phenotype data
pheno <- phen_dat

covs <- colnames(pheno)[!colnames(pheno) %in% c("aln", "alnqlet", "qlet", "Sample_Name", "ad")]

phen <- "ad"

# -------------------------------------------------------
# Generate SVs
# -------------------------------------------------------

svs <- generate_svs(trait = phen, 
					phen_data = pheno, 
					meth_data = mdata, 
					covariates = covs, 
					nsv = 10, 
					IID = "Sample_Name")

# -------------------------------------------------------
# Check association between SVs and phenotype of interest
# -------------------------------------------------------

sv_nam <- grep("sv", colnames(svs), value=T)

sv_check_dat <- pheno %>%
	left_join(svs) %>%
	dplyr::select(one_of(phen, sv_nam)) %>%
	na.omit

sv_assoc <- map_dfr(sv_nam, function(x) {
	print(x)
	form <- as.formula(paste(x, phen, sep = " ~ "))
	fit <- lm(form, data = sv_check_dat)
	out_nums <- summary(fit)$coef[2, ]
	out <- as_tibble(t(as.matrix(out_nums))) %>%
		mutate(sv = x) %>%
		dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`)
	return(out)
})

## remove associations at P<0.01 (changing from P<0.05 as doing 10 tests here...)
sv_to_rm <- sv_assoc %>% 
	dplyr::filter(P < 0.01) %>%
	pull(sv)

svs_out <- svs %>% 
	dplyr::select(-one_of(sv_to_rm))

# -------------------------------------------------------
# Write out results and info about removed SVs
# -------------------------------------------------------

write.table(svs_out, file = out_file,
			sep = "\t", quote = F, col.names = T, row.names = F)

sv_rem_info <- lapply(sv_to_rm, function(x) {svs[[x]]})
names(sv_rem_info) <- sv_to_rm
sv_rem_info$sv_assoc <- sv_assoc

save(sv_rem_info, file = removed_out_file)





