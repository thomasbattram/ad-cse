# ----------------------------------------
# Cell-specific EWAS of AD 
# ----------------------------------------

## Aim: Run EWAS using cell-specific effects methods for AD

## Date: 2022-08-12

## pkgs
library(tidyverse) # tidy code and data
library(matrixStats) # for imputing matrix
library(TCA) # TCA
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3]
covar_file <- args[4]
out_files <- args[5]
max_chunks <- as.numeric(args[6])

print(args)

# phen_file <- "../data-extraction-and-qc/data/ad-data-cleaned.tsv"
# meth_file <- "../data-extraction-and-qc/data/clean-meth.RData"
# svs_file <- ""
# covar_file <- "../data-extraction-and-qc/data/covars-no-cc.txt"
# out_files <- "results/ewas/tca-temp/tca-res-1.RData results/ewas/tca-temp/tca-res-2.RData"
# max_chunks <- 100

out_files <- unlist(str_split(out_files, " "))
# out_files <- paste0("results/ewas/tca-temp/tca-res-", 1:100, ".RData")

## read in data
pheno_dat <- read_tsv(phen_file)
meth <- new_load(meth_file)
# svs <- read_tsv(svs_file)

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

covs <- readLines(covar_file)
if (svs_file != "") {
    svs <- read_tsv(svs_file)
    pheno_dat <- pheno_dat %>%
        left_join(svs)       
}

phen_dat <- pheno_dat %>%
    dplyr::select(Sample_Name, aln, all_of(c(trait, covs))) %>%
    dplyr::filter(!is.na(Sample_Name)) %>%
    na.omit(.)

id_vars <- c("Sample_Name", "aln", "alnqlet", "qlet")
cell_types <- colnames(pheno_dat)[!colnames(pheno_dat) %in% c(id_vars, covs, trait)]
cell_counts <- as.matrix(pheno_dat[, c(cell_types)])
rownames(cell_counts) <- pheno_dat$Sample_Name

cell_counts <- cell_counts[rownames(cell_counts) %in% phen_dat$Sample_Name, ]
stopifnot(all(rownames(cell_counts) == phen_dat$Sample_Name))

meth <- meth[, phen_dat$Sample_Name]

## Make sex a 1/2 phenotype F=1, M=2
phen_dat$sex <- ifelse(phen_dat$sex == "F", 1, 2)

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

## make all categorical variables into binary variables
cat_covs <- c("BCD_plate")
for (cov in cat_covs) {
    if (!cov %in% colnames(phen_dat)) print("NEXT"); next
    uniq_vals <- unique(phen_dat[[cov]])
    for (i in seq_along(uniq_vals)) {
        val <- uniq_vals[i]
        phen_dat[[paste0(cov, i)]] <- ifelse(phen_dat[[cov]] == val, 1, 0)
    }
}
phen_dat <- dplyr::select(phen_dat, -one_of(cat_covs))

# ----------------------------------------
# EWAS functions
# ----------------------------------------

run_tca <- function(temp_phen, temp_meth, phen, cc, covs, IID)
{
    phen_vals <- temp_phen[[phen]]
    tca_phen <- matrix(phen_vals)
    colnames(tca_phen) <- phen
    rownames(tca_phen) <- temp_phen[[IID]]
    tca_covs <- as.matrix(temp_phen[, !colnames(temp_phen) %in% c(phen, IID, "aln")])
    rownames(tca_covs) <- temp_phen[[IID]]
    out <- tca(X = temp_meth, C1 = tca_phen, W = cc, C2 = tca_covs)
    return(out)
}

sort_tca <- function(tca_res) 
{
    tca_beta <- tca_res$gammas_hat
    tca_p <- tca_res$gammas_hat_pvals
    tca_out <- list(beta = tca_beta[, grep(trait, colnames(tca_beta))], 
                      p = tca_p[, grep(trait, colnames(tca_p))])
    tca_out <- lapply(tca_out, function(x) {
        colnames(x) <- gsub("\\..*", "", colnames(x))
        return(x)
    })
    return(tca_out)
}

run_ewas <- function(phen, p_dat, cc, meth_dat, IID, method, covs) 
{
    # Match meth to Pheno
    meth_dat <- meth_dat[, colnames(meth_dat) %in% p_dat[[IID]]]
    meth_dat <- meth_dat[, na.omit(match(p_dat[[IID]], colnames(meth_dat)))]
    temp_phen <- p_dat[match(colnames(meth_dat), p_dat[[IID]]), ]
    temp_cc <- cc[rownames(cc) %in% temp_phen$Sample_Name, ]
    
    if (!all(temp_phen[[IID]] == colnames(meth_dat))) stop("phenotype and DNAm data not matched.")
    ## FOR TESTS
    # temp_meth <- temp_meth[1:50, 1:50]
    # temp_phen <- temp_phen[1:50, ]
    # temp_cc <- temp_cc[1:50, ]

    function_name <- paste0("run_", method)
    ewas_func <- match.fun(function_name)
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

## test
# meth <- meth[1:10000,]

## subset meth matrix 
chunk_lengths <- round(seq(1, nrow(meth), length = max_chunks+1))
if (!nrow(meth) %in% chunk_lengths) stop("Not going across all CpG sites")

method <- "tca"
# out_file <- out_files[1]
lapply(out_files, function(out_file) {
    chunk <- as.numeric(str_extract(out_file, "[0-9][0-9]*"))
    meth_dat <- meth[chunk_lengths[chunk]:chunk_lengths[chunk+1], ]
    message("Running TCA analyses using chunk: ", chunk)
    ewas_res <- run_ewas(phen = trait, 
                         p_dat = phen_dat, 
                         cc = cell_counts, 
                         meth_dat = meth_dat,
                         IID = "Sample_Name", 
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
    rm(list = c("out_res", "ewas_res"))    
})

print("FIN")


