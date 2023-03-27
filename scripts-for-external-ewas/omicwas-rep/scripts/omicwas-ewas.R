# ----------------------------------------
# Cell-specific EWAS of AD using omicWAS 
# ----------------------------------------

## Aim: Run EWAS using omicWAS for AD

## pkgs
library(tidyverse) # tidy code and data
library(omicwas) # for running omicwas EWAS
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3] # make this "" if not using
covar_file <- args[4]
hits_file <- args[5]
out_file <- args[6]

# phen_file <- "data/cleaned-pheno-data.tsv"
# meth_file <- "data/clean-meth.RData"
# svs_file <- "data/svs/ad-svs.tsv"
# covar_file <- "data/covars-no-cc.txt"
# hits_file <- "data/meta-analysis-hits.txt"
# out_file <- "results/ewas/omicwas-res.RData"

message("Here are your args: \n", 
        "phen_file = ", phen_file, "\n", 
        "meth_file = ", meth_file, "\n", 
        "svs_file = ", svs_file, "\n", 
        "covar_file = ", covar_file, "\n", 
        "hits_file = ", hits_file, "\n", 
        "out_file = ", out_file, "\n", 
        )

## read in data
pheno_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)
prev_hits <- readLines(hits_file)
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

## Get only CpGs of interest
meth <- meth[rownames(meth) %in% prev_hits, ]

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
 
run_omicwas <- function(temp_phen, temp_meth, phen, cc, covs, IID, seed = 2)
{
    phen_vals <- temp_phen[[phen]]
    omicwas_phen <- matrix(phen_vals)
    colnames(omicwas_phen) <- phen
    rownames(omicwas_phen) <- temp_phen[[IID]]
    omicwas_covs <- as.matrix(temp_phen[, !colnames(temp_phen) %in% c(phen, IID, "aln")])
    rownames(omicwas_covs) <- temp_phen[[IID]]
    res <- ctassoc(X = omicwas_phen, W = cc, Y = temp_meth, C = omicwas_covs, 
                   test = "nls.logit", regularize = TRUE)
    ## Matrix (or vector) of covariates; samples x covariates. X, W, Y, C should benumeric.
    out <- res$coefficients %>% dplyr::filter(term == phen)
    return(out) 
}

sort_omicwas <- function(omicwas_res, prev_hits)
{
    cols_to_get <- c("estimate", "p.value")
    omicwas_out <- lapply(cols_to_get, function(x) {
        out <- omicwas_res %>%
            dplyr::select(response, celltype, one_of(x)) %>%
            pivot_wider(names_from = celltype, values_from = one_of(x)) %>%
            as.data.frame
        rownames(out) <- out$response
        out <- out[, !colnames(out) == "response"]
        return(out)
    })
    names(omicwas_out) <- c("beta", "p")
    return(omicwas_out)
}

run_ewas <- function(phen, p_dat, cc, meth_dat, IID, method, covs) 
{
    # Match meth to Pheno
    temp_meth <- meth_dat[, colnames(meth_dat) %in% p_dat[[IID]]]
    temp_meth <- meth_dat[, na.omit(match(p_dat[[IID]], colnames(meth_dat)))]
    temp_phen <- p_dat[match(colnames(temp_meth), p_dat[[IID]]), ]
    temp_cc <- cc[rownames(cc) %in% temp_phen$Sample_Name, ]
    
    if (!all(temp_phen[[IID]] == colnames(temp_meth))) stop("phenotype and DNAm data not matched.")
    ## FOR TESTS
    # temp_meth <- temp_meth[1:50, 1:50]
    # temp_phen <- temp_phen[1:50, ]
    # temp_cc <- temp_cc[1:50, ]

    function_name <- paste0("run_", method)
    ewas_func <- match.fun(function_name)
    # p_to_keep <- p_to_keep[1]
    message("EWAS TIME")
    res <- tryCatch({
        ewas_func(temp_phen, temp_meth, phen, temp_cc, covs, IID)
    }, error = function(e) {
        usr_m <- paste0("Error in EWAS of ", phen, " using ", method, ".")
        err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = phen)        
    })
    return(res)
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------

method <- "omicwas"
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
out_res <- sort_omicwas(ewas_res)
message("Saving sorted results to ", out_file)
save(out_res, file = out_file)

print("FIN")
