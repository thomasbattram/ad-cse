# -------------------------------------------------------
# Script to generate surrogate variables
# -------------------------------------------------------

## pkgs
library(tidyverse) # tidy data and code
library(sva) # calculating surrogate variables
library(SmartSVA) # calculating SVs
library(usefunc) # useful utility functions

## args
args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
batch_vars_file <- args[2]
meth_file <- args[3]
out_file <- args[4]
removed_out_file <- args[5]
covars_outfile <- args[6]
cc_covars_outfile <- args[7]
heatmap_outfile <- args[8]

## Manual args
# phen_file <- "data/cleaned-pheno-data.tsv"
# batch_vars_file <- "data/batch-vars.txt"
# meth_file <- "data/clean-meth.RData"
# out_file <- "data/svs/ad-svs.tsv"
# removed_out_file <- "data/svs/ad-removed-svs.RData"
# covars_outfile <- "data/covars-no-cc.txt"
# cc_covars_outfile <- "data/covars-cc.txt"
# heatmap_outfile <- "results/svs-heatmap.png"

message("Here are your args: \n", 
		"phen_file = ", phen_file, "\n", 
		"batch_vars_file = ", batch_vars_file, "\n",  
		"meth_file = ", meth_file, "\n", 
		"out_file = ", out_file, "\n", 
		"removed_out_file = ", removed_out_file, "\n", 
		"covars_outfile = ", covars_outfile, "\n", 
		"cc_covars_outfile = ", cc_covars_outfile, "\n", 
		"heatmap_outfile = ", heatmap_outfile, "\n", 
		)

source(useful_funcs_script)

## read in data
pheno <- read_tsv(phen_file)
meth <- new_load(meth_file)
batch_vars <- readLines(batch_vars_file)
cc_covars <- readLines(cc_covars_outfile)
no_cc_covars <- readLines(covars_outfile)

# -------------------------------------------------------
# Check data
# -------------------------------------------------------

expected_cols <- c("IID", "ad", batch_vars, cc_covars)
lapply(expected_cols, function(col) {
	if (!col %in% colnames(phen_dat)) {
		stop(col, " not in the column names of your phenotype data. ",
				  "Please make sure all batch variables, the identifier variable (IID), the covariates, ", 
				  "and the AD case-control variable (ad) are within the data.frame 'pheno'")
	}
})

# -------------------------------------------------------
# functions
# -------------------------------------------------------

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
generate_svs <- function(trait, phen_data, mdat, covariates = "", nsv, 
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
	
	mdat <- mdat[, colnames(mdat) %in% phen[[IID]]]
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

	message("mdat dimensions: ", paste0(nrow(mdat), ", ", ncol(mdat)))
	message("pheno dimensions: ", paste0(nrow(phen), ", ", ncol(phen)))

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

# phenotype data
phen <- "ad"
covs <- cc_covars

## methylation data
mdata <- impute_matrix(meth)
rm(meth)

## adjust methylation data for cell counts and take residuals
celltypes <- covs[!covs %in% no_cc_covars]
form <- paste(paste0("pheno$", celltypes), collapse = " + ")
ori_time <- proc.time()
meth_resid <- lapply(1:nrow(mdata), function(x) {
	print(x)
	dnam <- mdata[x,, drop=T]
	resid(lm(as.formula(paste0("dnam ~ ", form))))
})
names(meth_resid) <- rownames(mdata)
mdata <- do.call(rbind, meth_resid)
message("mdata dimensions: ", paste0(nrow(mdata), ", ", ncol(mdata)))
time_taken <- proc.time() - ori_time
time_taken # ~ 41 mins
rm(meth_resid)

# -------------------------------------------------------
# Generate SVs
# -------------------------------------------------------

svs <- generate_svs(trait = phen, 
					phen_data = pheno, 
					mdat = mdata, 
					covariates = covs[!covs %in% celltypes], 
					nsv = 10, 
					IID = "IID")

# svs <- read_tsv(out_file)

# -------------------------------------------------------
# Check association between SVs and phenotype of interest
# -------------------------------------------------------

sv_nam <- grep("sv", colnames(svs), value=T)

#' Assess association with SVs
#' 
#' @param svs table of SV values
#' @param pheno phenotype data
#' @param trait trait of interest (cov or phenotype)
regress_sv <- function(svs, pheno, trait)
{
	sv_check_dat <- pheno %>%
		left_join(svs) %>%
		dplyr::select(one_of(trait, sv_nam)) %>%
		na.omit

	sv_assoc <- map_dfr(sv_nam, function(x) {
		print(x)
		form <- as.formula(paste(x, trait, sep = " ~ "))
		fit <- lm(form, data = sv_check_dat)
		out_nums <- summary(fit)$coef[2, ]
		out_r2 <- summary(fit)$adj.r.squared
		out <- as_tibble(t(as.matrix(out_nums))) %>%
			mutate(sv = x, r2 = out_r2) %>%
			dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`, r2)
		return(out)
	})

	## remove associations at P<0.01 (changing from P<0.05 as doing 10 tests here...)
	sv_to_rm <- sv_assoc %>% 
		dplyr::filter(P < 0.01) %>%
		pull(sv)

	out <- list(assoc = sv_assoc, to_rm = sv_to_rm)
	return(out)
}

sv_vars <- c(phen, covs, batch_vars)
sv_assoc_list <- lapply(sv_vars, regress_sv, svs = svs, pheno = pheno)
names(sv_assoc_list) <- sv_vars

sv_assoc <- bind_rows(map(sv_assoc_list, "assoc"), .id = "variable")

## reorder the variables for the plot
sv_assoc$sv <- factor(sv_assoc$sv, levels = sv_nam)
var_levels <- c(sv_vars[!sv_vars %in% celltypes], celltypes)
sv_assoc$variable <- factor(sv_assoc$variable, levels = var_levels)

## Table
# SV | trait | beta

heatmap_r2 <- ggplot(sv_assoc, aes(sv, variable, fill = r2)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = round(r2, 2))) + 
    scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="r2")

ggsave(heatmap_outfile, plot = heatmap_r2)

## REMOVE SVS ASSOC WITH TRAIT HERE
sv_to_rm <- sv_assoc %>%
	dplyr::filter(variable == "ad") %>%
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
sv_rem_info$sv_assoc <- sv_assoc %>%
	dplyr::filter(variable == phen)

save(sv_rem_info, file = removed_out_file)

## Add SVs to covariate file
sv_covs <- grep("sv", colnames(svs_out), value=T)
message("Appending the following variables to the covariates file: ", paste(sv_covs, collapse = ", "))
covars <- readLines(covars_outfile)
covars <- unique(c(covars, sv_covs))
write(covars, file = covars_outfile)
covars_cc <- readLines(cc_covars_outfile)
covars_cc <- unique(c(covars_cc, sv_covs))
write(covars_cc, file = cc_covars_outfile)



