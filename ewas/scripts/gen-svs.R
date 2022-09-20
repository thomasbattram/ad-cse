# -------------------------------------------------------
# Script to generate surrogate variables for ARIES
# -------------------------------------------------------

# srun --job-name "InteractiveJob" --partition=veryshort --nodes=1 --ntasks-per-node=4 --cpus-per-task=4 --time=6:00:00 --mem=50GB --pty bash

## pkgs
library(tidyverse) # tidy data and code
library(sva) # calculating surrogate variables
library(aries) # get aries data - for batch variables
library(SmartSVA) # calculating SVs
library(matrixStats) # imputing DNAm data
library(usefunc) # personal package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
aries_dir <- args[3]
out_file <- args[4]
removed_out_file <- args[5]
covars_outfile <- args[6]
cc_covars_outfile <- args[7]
heatmap_outfile <- args[8]

# phen_file <- "../data-extraction-and-qc/data/ad-data-cleaned.tsv"
# meth_file <- "../data-extraction-and-qc/data/clean-meth.RData"
# aries_dir <- "/user/work/ms13525/aries"
# out_file <- "data/svs/ad-svs.tsv"
# removed_out_file <- "data/svs/ad-removed-svs.RData"
# covars_outfile <- "../data-extraction-and-qc/data/covars-no-cc.txt"
# cc_covars_outfile <- "../data-extraction-and-qc/data/covars-cc.txt"
# heatmap_outfile <- "results/svs-heatmap.png"

## read in data
phen_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)
aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
batch_vars <- c("BCD_plate", "Slide")
phen_dat <- phen_dat %>%
	left_join(aries$samples[, c("Sample_Name", batch_vars)])
rm(aries)

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
pheno <- phen_dat
rm(phen_dat)
covs <- colnames(pheno)[!colnames(pheno) %in% c("aln", "alnqlet", "qlet", "Sample_Name", "ad", batch_vars)]

phen <- "ad"

## methylation data
mdata <- impute_matrix(meth)
rm(meth)

## adjust methylation data for cell counts and take residuals
## FOR TESTING
# mdata <- mdata[sample(1:nrow(mdata), 5e4), ]
###
celltypes <- covs[!covs %in% c("age", "sex")]
form <- paste(paste0("pheno$", celltypes), collapse = " + ")
# test_mdata <- mdata[1:100,]
ori_time <- proc.time()
meth_resid <- lapply(1:nrow(mdata), function(x) {
	print(x)
	dnam <- mdata[x,, drop=T]
	resid(lm(as.formula(paste0("dnam ~ ", form))))
})
# names(meth_resid) <- rownames(test_mdata)
# test_mdata <- do.call(rbind, meth_resid)
names(meth_resid) <- rownames(mdata)
mdata <- do.call(rbind, meth_resid)
message("mdata dimensions: ", paste0(nrow(mdata), ", ", ncol(mdata)))
# mdata <- t(mdata)
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
					IID = "Sample_Name")

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
		out <- as_tibble(t(as.matrix(out_nums))) %>%
			mutate(sv = x) %>%
			dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`)
		return(out)
	})

	## remove associations at P<0.01 (changing from P<0.05 as doing 10 tests here...)
	sv_to_rm <- sv_assoc %>% 
		dplyr::filter(P < 0.01) %>%
		pull(sv)

	out <- list(assoc = sv_assoc, to_rm = sv_to_rm)
	return(out)
}

sv_vars <- colnames(pheno)[!colnames(pheno) %in% c("aln", "alnqlet", "qlet", "Sample_Name")]
sv_assoc_list <- lapply(sv_vars, regress_sv, svs = svs, pheno = pheno)
names(sv_assoc_list) <- sv_vars

sv_assoc <- bind_rows(map(sv_assoc_list, "assoc"), .id = "variable")

## Table
# SV | trait | beta

heatmap_b <- ggplot(sv_assoc, aes(sv, variable, fill = beta)) +
    geom_tile(color = "white") + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, space = "Lab", 
                         name="Beta")

ggsave(heatmap_outfile, plot = heatmap_b)

## START FROM HERE!!

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
	dplyr::filter(variable == "ad")

save(sv_rem_info, file = removed_out_file)

## Add SVs to covariate file
sv_covs <- grep("sv", colnames(svs_out), value=T)
message("Appending the following variables to the covariates file: ", paste(sv_covs, collapse = ", "))
write(sv_covs, file = covars_outfile, append = TRUE)
write(sv_covs, file = cc_covars_outfile, append = TRUE)



