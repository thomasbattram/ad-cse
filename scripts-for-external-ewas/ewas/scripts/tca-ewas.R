# ----------------------------------------
# Cell-specific EWAS of AD using TCA
# ----------------------------------------

## Aim: Run EWAS using cell-specific effects methods for AD

## pkgs
library(tidyverse) # tidy code and data
library(TCA) # TCA
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
svs_file <- args[3]
covar_file <- args[4]
out_files <- args[5]
max_chunks <- as.numeric(args[6])

## manual args
# phen_file <- "data/cleaned-pheno-data.tsv"
# meth_file <- "data/clean-meth.RData"
# svs_file <- "data/svs/ad-svs.tsv"
# covar_file <- "data/covars-no-cc.txt"
# out_files <- paste0(paste0("results/ewas/tca-temp/tca-res-", 1:100, ".RData"), collapse = " ")
# max_chunks <- 100

message("Here are your args: \n", 
        "phen_file = ", phen_file, "\n", 
        "meth_file = ", meth_file, "\n", 
        "svs_file = ", svs_file, "\n", 
        "covar_file = ", covar_file, "\n", 
        "out_files = ", out_files, "\n", 
        "max_chunks = ", max_chunks
        )


out_files <- unlist(str_split(out_files, " "))

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
    filename <- basename(out_file)
    chunk <- as.numeric(str_extract(filename, "[0-9][0-9]*"))
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


