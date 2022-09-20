# ----------------------------------------
# ewaff EWAS script
# ----------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # for EWAS functions
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3] # make this "" if not using
covar_file <- args[4]
out_file <- args[5]

# phen_file <- "../data-extraction-and-qc/data/ad-data-cleaned.tsv"
# meth_file <- "../data-extraction-and-qc/data/clean-meth.RData"
# svs_file <- ""
# covar_file <- "../data-extraction-and-qc/data/covars-no-cc.txt"
# out_file <- "results/ewas/ewaff-res-no-cc.tsv"

model <- ifelse(grepl("no-cc", out_file), "no_cc", "cc")

## read in data
phen_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)

options(mc.cores=1)

# ----------------------------------------
# EWAS functions
# ----------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, IID, covs)
{
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

    return(temp_phen)
}

run_ewas <- function(phen, pheno_dat, svs_file, meth_dat, IID, out_file, covs) 
{
    # prep pheno data
    temp_phen <- prep_pheno_data(phen, pheno_dat, svs_file, IID, covs)
    sv_nam <- grep("sv", colnames(temp_phen), value = T)
    covs <- c(covs, sv_nam)
    covs <- covs[covs %in% colnames(temp_phen)]

    # Match meth to Pheno
    temp_meth <- meth_dat[, na.omit(match(temp_phen[[IID]], colnames(meth_dat)))]
    temp_phen <- temp_phen[match(colnames(temp_meth), temp_phen[[IID]]), ]

    # Get cases and controls per probe
    cases <- temp_phen[temp_phen[[phen]] == 1, IID, drop = TRUE]

	n_cases <- rowSums(!is.na(temp_meth[, cases]))
	n_controls <- rowSums(!is.na(temp_meth[, !colnames(temp_meth) %in% cases]))
	probe_n_cc <- as_tibble(cbind(probeID = names(n_cases), N_cases = n_cases, N_controls = n_controls))

    if (!all(temp_phen[[IID]] == colnames(temp_meth))) stop("phenotype and DNAm data not matched.")

    model <- as.formula(paste0(phen, " ~ ", paste(c("methylation", covs), collapse = " + ")))

    print(temp_phen)

    # Run EWAS using ewaff
    message("Running EWAS")
    obj <- tryCatch({
        obj <- ewaff.sites(model, variable.of.interest = phen,
                           methylation = temp_meth, data = temp_phen, method = "glm", 
                           generate.confounders = NULL, family = "binomial")
    }, error = function(e) {
        usr_m <- paste0("Error in EWAS of ", phen)
        err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = phen)
    })
    # free up some space
    rm(temp_meth)

    if (length(obj) == 1) {
        return(NULL)
    }
    res <- obj$table %>%
        rownames_to_column(var = "probeID") %>%
        dplyr::select(probeID, BETA = estimate, SE = se, P = p.value, N = n) %>%
        left_join(probe_n_cc, by = "probeID")

    write.table(res, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F)
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------
phen <- "ad"

covariates <- readLines(covar_file)

run_ewas(phen = phen, 
		 pheno_dat = phen_dat, 
		 svs_file = svs_file,
		 meth_dat = meth,
		 IID = "Sample_Name", 
		 out_file = out_file, 
		 covs = covariates)

print("FIN")