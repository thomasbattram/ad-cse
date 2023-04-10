# ----------------------------------------
# Joint location-and-scale score test EWAS script
# ----------------------------------------

## This script uses the "jlst" R package to run join location-and-scale score test
## to see if we can detect mean or variance differences in DNAm between AD cases and controls 
## when jointly modelling them.

## pkgs
library(tidyverse) # tidy code and data
## devtools::install_github("jrs95/jlst")
library(jlst) # functions for variance tests
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3] # make this "" if not using
covar_file <- args[4]
out_file <- args[5]

message("Arguments are: ", args)

# phen_file <- "../data-extraction-and-qc/data/ad-data-cleaned.tsv"
# meth_file <- "../data-extraction-and-qc/data/clean-meth.RData"
# svs_file <- "data/svs/ad-svs.tsv"
# covar_file <- "../data-extraction-and-qc/data/covars-cc.txt"
# out_file <- "results/ewas/jlssc-res-cc.tsv"

use_covs <- ifelse(grepl("cc", out_file), TRUE, FALSE)

message("It is ", use_covs, " that covariates will be used in this analysis")

## read in data
phen_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)
covariates <- readLines(covar_file)

# meth <- meth[1:100, ]

# ----------------------------------------------------------------
# ewas functions
# ----------------------------------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, IID, covs)
{
	message("Preparing phenotype data")
    if (svs_file == "") return(na.omit(pheno_dat[, c(IID, phen, covs)]))
    ## read in svs and cell counts
    svs <- read_tsv(svs_file)
    sv_nam <- grep("sv", colnames(svs), value = T)
    covs <- c(covs, sv_nam)
    # Prepare phenotype data
    temp_phen <- pheno_dat %>%
        left_join(svs, by = setNames(IID, IID)) %>%
        dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
        na.omit(.)

    message("Phenotype data prepared")

    return(temp_phen)
}

varewas <- function(use_covs, cpg, temp_phen, temp_meth, covs) {
	message("Running vartest")
	covar_dat <- as.data.frame(temp_phen[, covs])
	if(!use_covs) covar_dat <- NULL
	var_out <- tryCatch({
		jlssc(y = temp_meth[rownames(temp_meth) == cpg, , drop = TRUE], 
			  x = as.factor(temp_phen[[phen]]), 
			  covar = covar_dat, 
			  type = 2)
    }, error = function(e) {
        usr_m <- paste0("Error in variance test of ", cpg)
        err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = cpg)
    })
    return(var_out)
}

calc_increments <- function(length_task, percent_increments = 10) 
{
	### Calculate increments to output message
	increments <- round(length_task / percent_increments)
	message_at <- round(seq(increments, length_task, length.out = percent_increments))
	names(message_at) <- seq(percent_increments, 100, percent_increments)
	return(message_at)
}

output_percent_complete <- function(n_task, increments) 
{
	### output percentage complete
	if (n_task %in% increments) {
		percent <- increments[increments == n_task]
		message(names(percent), "% complete.")
	}
}


run_ewas <- function(phen, pheno_dat, svs_file, meth_dat, IID, out_file, covs) 
{
    # prep pheno data
    temp_phen <- prep_pheno_data(phen, pheno_dat, svs_file, IID, covs)
    # sv_nam <- grep("sv", colnames(temp_phen), value = T)
    # covs <- c(covs, sv_nam)
    covs <- covs[covs %in% colnames(temp_phen)]

    # Match meth to Pheno
    temp_meth <- meth_dat[, na.omit(match(temp_phen[[IID]], colnames(meth_dat)))]
    temp_phen <- temp_phen[match(colnames(temp_meth), temp_phen[[IID]]), ]

    # Get cases and controls per probe
    cases <- temp_phen[temp_phen[[phen]] == 1, IID, drop = TRUE]

	n_cases <- rowSums(!is.na(temp_meth[, cases]))
	n_controls <- rowSums(!is.na(temp_meth[, !colnames(temp_meth) %in% cases]))
	probe_n_cc <- as_tibble(cbind(probeID = names(n_cases), N = n_cases + n_controls, N_cases = n_cases, N_controls = n_controls))

    if (!all(temp_phen[[IID]] == colnames(temp_meth))) stop("phenotype and DNAm data not matched.")

    print(temp_phen)

	increments <- calc_increments(length_task = nrow(temp_meth))

    # Run variance EWAS using jlst
    message("Running EWAS")
    list_res <- lapply(1:nrow(temp_meth), function(x) {
    	output_percent_complete(x, increments)
    	cpg <- rownames(temp_meth)[x]
    	res <- varewas(use_covs, cpg, temp_phen, temp_meth, covs)
    	if (class(res) == "character") return(res)
    	res$probeID <- cpg
    	return(res)
    })
    
    # free up some space
    rm(temp_meth)

    succ_res <- list_res[sapply(list_res, class) != "character"]

    res <- bind_rows(succ_res) %>%
    	left_join(probe_n_cc, by = "probeID") %>%
    	dplyr::select(probeID, everything())
    write.table(res, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F)
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------
phen <- "ad"

run_ewas(phen = phen, 
		 pheno_dat = phen_dat, 
		 svs_file = svs_file,
		 meth_dat = meth,
		 IID = "Sample_Name", 
		 out_file = out_file, 
		 covs = covariates)

print("FIN")