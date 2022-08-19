# ----------------------------------------
# Combine TCA EWAS results
# ----------------------------------------

## Aim: Combine res from TCA analyses

## Date: 2022-08-19

## pkgs
library(tidyverse) # tidy code and data
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
in_files <- args[1]
out_file <- args[2]

# in_files <- "results/ewas/tca-temp/tca-res-1.RData results/ewas/tca-temp/tca-res-2.RData"
# out_file <- "results/ewas/tca-res.RData"

in_files <- unlist(str_split(in_files, " "))
all_res <- lapply(in_files, function(file) {out <- new_load(file); return(out)})

out_res_names <- c("beta", "p")
out_res <- lapply(out_res_names, function(x) {
	res <- map(all_res, x)
	out <- do.call(rbind, res)
	out <- out[!duplicated(rownames(out)), ]
	if (any(duplicated(rownames(out)))) stop("Duplicated rownames in TCA combined output")
	return(out)
})
names(out_res) <- out_res_names
str(out_res)

save(out_res, file = out_file)