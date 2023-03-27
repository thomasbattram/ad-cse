# ----------------------------------------
# Cell-specific EWAS of AD 
# ----------------------------------------

## Aim: Run EWAS using cell-specific effects methods for AD

## pkgs
library(tidyverse) # tidy code and data
library(matrixStats) # for imputing matrix
library(EpiDISH) # celldmc
library(usefunc)

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3] # make this "" if not using
covar_file <- args[4]
out_file <- args[5]

## manual args
# phen_file <- "data/cleaned-pheno-data.tsv"
# meth_file <- "data/clean-meth.RData"
# svs_file <- "data/svs/ad-svs.tsv"
# covar_file <- "data/covars-no-cc.txt"
# out_file <- "results/ewas/celldmc-res.RData"

message("Here are your args: \n", 
        "phen_file = ", phen_file, "\n", 
        "meth_file = ", meth_file, "\n", 
        "svs_file = ", svs_file, "\n", 
        "covar_file = ", covar_file, "\n", 
        "out_file = ", out_file, "\n", 
        )

## read in data
pheno_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)
covs <- readLines(covar_file)
if (svs_file != "") {
    svs <- read_tsv(svs_file)
    pheno_dat <- pheno_dat %>%
        left_join(svs)       
}
trait <- "ad"
IID <- "IID"

# -------------------------------------------------------
# sort data for running EWAS
# -------------------------------------------------------

meth <- impute_matrix(meth)

phen_dat <- pheno_dat %>%
    dplyr::select(all_of(c(IID, trait, covs))) %>%
    dplyr::filter(!is.na(IID)) %>%
    na.omit(.)

cell_types <- colnames(pheno_dat)[!colnames(pheno_dat) %in% c(IID, covs, trait)]
cell_counts <- as.matrix(pheno_dat[, c(cell_types)])
rownames(cell_counts) <- pheno_dat[[IID]]

cell_counts <- cell_counts[rownames(cell_counts) %in% phen_dat[[IID]], ]
stopifnot(all(rownames(cell_counts) == phen_dat[[IID]]))

meth <- meth[, phen_dat[[IID]]]

## Some methods can't be used with negative cell count estimates. These should be very close to 0 anyway, so setting to 0
cell_counts[sign(cell_counts) == -1] <- 0

## Need to re-estimate cell props so they all add up to 1 for TCA to work so doing it for all methods
reest_cell_props <- function(cc_mat)
{
    out_mat <- map_dfr(1:nrow(cc_mat), function(x) {
        cc_mat[x, ] / sum(cc_mat[x, ])
    })
    return(as.matrix(out_mat))
}

rownam <- rownames(cell_counts)
cell_counts <- reest_cell_props(cell_counts)
rownames(cell_counts) <- rownam

# ----------------------------------------
# EWAS functions
# ----------------------------------------

run_celldmc <- function(temp_phen, temp_meth, phen, cc, covs, IID)
{
    cov_form <- as.formula(paste0("~ ", paste(covs, collapse = " + ")))
    cov_mat <- model.matrix(cov_form, data = temp_phen)
    out <- CellDMC(beta.m = temp_meth, pheno.v = temp_phen[[phen]], frac.m = cc, cov.mod = cov_mat)
    return(out)
}

sort_celldmc <- function(celldmc_res) 
{
    celldmc_coefs <- celldmc_res$coe
    cols_to_get <- c("Estimate", "SE", "p")
    celldmc_out <- lapply(cols_to_get, function(x) {
        out <- map_dfr(celldmc_coefs, x) %>%
            as.data.frame
        rownames(out) <- rownames(celldmc_coefs[[1]])
        return(out)
    })
    names(celldmc_out) <- c("beta", "se", "p")
    return(celldmc_out)
}
 
run_ewas <- function(phen, p_dat, cc, meth_dat, IID, method, covs) 
{
    # Match meth to Pheno
    meth_dat <- meth_dat[, colnames(meth_dat) %in% p_dat[[IID]]]
    meth_dat <- meth_dat[, na.omit(match(p_dat[[IID]], colnames(meth_dat)))]
    temp_phen <- p_dat[match(colnames(meth_dat), p_dat[[IID]]), ]
    temp_cc <- cc[rownames(cc) %in% temp_phen[[IID]], ]
    
    if (!all(temp_phen[[IID]] == colnames(meth_dat))) stop("phenotype and DNAm data not matched.")
    ## FOR TESTS
    # meth_dat <- meth_dat[1:50, 1:50]
    # temp_phen <- temp_phen[1:50, ]
    # temp_cc <- temp_cc[1:50, ]

    function_name <- paste0("run_", method)
    ewas_func <- match.fun(function_name)
    # p_to_keep <- p_to_keep[1]

    # print(colSums(model.matrix(~ temp_cc + temp_phen[[phen]]:temp_cc)[, -1]) == 0)

    message("EWAS TIME")
    res <- tryCatch({
        ewas_func(temp_phen, meth_dat, phen, temp_cc, covs, IID)
    }, error = function(e) {
        usr_m <- paste0("Error in EWAS of ", phen, " using ", method, ".")
        err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = phen)        
    })
    return(res)
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------

method <- "celldmc"
message("Running analyses using ", method)
ewas_res <- run_ewas(phen = trait, 
                     p_dat = phen_dat, 
                     cc = cell_counts, 
                     meth_dat = meth,
                     IID = IID, 
                     method = method, 
                     covs = covs)
if (is.character(ewas_res)) {
    stop("EWAS of AD using ", method, " failed to run.")
}
message("Finished analyses using ", method)
sort_function_name <- paste0("sort_", method)
sort_func <- match.fun(sort_function_name)
out_res <- sort_func(ewas_res)
message("Saving sorted results to ", out_file)
save(out_res, file = out_file)

print("FIN")
