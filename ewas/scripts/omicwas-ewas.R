# ----------------------------------------
# Cell-specific EWAS of AD 
# ----------------------------------------

## Aim: Run EWAS using cell-specific effects methods for AD

## Date: 2022-08-12

## pkgs
library(tidyverse) # tidy code and data
library(matrixStats) # for imputing matrix
library(omicwas) # for running omicwas EWAS
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3]
hits_file <- args[4]
out_file <- args[5]

# phen_file <- "../data-extraction-and-qc/data/ad-data-cleaned.tsv"
# meth_file <- "../data-extraction-and-qc/data/clean-meth.RData"
# svs_file <- "data/svs/ad-svs.tsv"
# hits_file <- "results/celldmc-tca-hits.RData"
# out_file <- "results/ewas/omicwas-res.RData"

## read in data
pheno_dat <- read_tsv(phen_file)
prev_hits <- new_load(hits_file)
meth <- new_load(meth_file)
svs <- read_tsv(svs_file)

trait <- "ad"

# -------------------------------------------------------
# sort data for running EWAS
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

meth <- impute_matrix(meth)

## Get only CpGs of interest
prev_hits <- prev_hits$res
cpgs <- unlist(map(prev_hits, "CpG"))
meth <- meth[rownames(meth) %in% cpgs, ]

cell_types <- c("Bcell", "CD4T", "CD8T", "Eos", "Mono", "Neu", "NK") ## ADD TO ME 
cell_counts <- as.matrix(pheno_dat[, c(cell_types)])
rownames(cell_counts) <- pheno_dat$Sample_Name

covs <- c(grep("sv", colnames(svs), value = T)) ## ADD MORE TO ME PLEASE
phen_dat <- pheno_dat %>%
    left_join(svs) %>%
    dplyr::select(Sample_Name, aln, all_of(c(trait, covs))) %>%
    dplyr::filter(!is.na(Sample_Name)) %>%
    na.omit(.)

cell_counts <- cell_counts[rownames(cell_counts) %in% phen_dat$Sample_Name, ]
stopifnot(all(rownames(cell_counts) == phen_dat$Sample_Name))

## Negative values are soooo close to zero, and values need to be 0-1 for some methods
cell_counts[sign(cell_counts) == -1] <- 0

## Need to re-estimate cell props so they all add up to 1 for TCA to work...
reest_cell_props <- function(cc_mat)
{
    out_mat <- map_dfr(1:nrow(cc_mat), function(x) {
        cc_mat[x,]/sum(cc_mat[x,])
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
                     IID = "Sample_Name", 
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
