# ----------------------------------------------------------------------
# Get f24 measured directly measured cell counts
# ----------------------------------------------------------------------

## pkgs
library(alspac) # extract data
library(tidyverse) # tidy code and data
library(openxlsx) # write out data as an excel spreadsheet
library(sftp)
library(labelled) # label data
library(usefunc) # personal package of useful functions

outfile <- "f24-direct-cellcounts.tsv"
outdir <- "."

# ----------------------------------------------------------------------
# sort alspac
# ----------------------------------------------------------------------

alspac_data_dir <- "/Volumes/Data/"
setDataDir(alspac_data_dir)

data(current)

# data(useful)

## RUN THIS TO CREATE THE DICTIONARIES
# current <- createDictionary("Current", name="current")
# useful <- createDictionary("Useful_data", name="useful")

## RUN THIS TO UPDATE THE DICTIONARIES
# updateDictionaries()

# ----------------------------------------------------------------------
# Extract variables and labels
# ----------------------------------------------------------------------

vars <- c("Neutrophils_F24", "Lymphocytes_F24", "Monocytes_F24", "Eosinophils_F24",
		  "Basophils_F24")

new_current <- current %>%
	dplyr::filter(name %in% vars)

## extraction of data
result <- extractVars(new_current)

## extraction of labels
val_labs <- usefunc::extract_alspac_labels(new_current, alsp_dir = alspac_data_dir)
unique_labs <- unique(unlist(val_labs))

## CHECK DATA HERE!!
test <- result[1:10, 1:10]
new_res <- result
for (ul in unique_labs) {
	new_res[new_res == ul] <- NA 
}

new_res <- new_res %>%
	dplyr::select(aln, qlet, alnqlet, all_of(vars))

summary(new_res)

## There is still -10 in there, that must be an artifact so removing
new_res[new_res == -10] <- NA

summary(new_res)

# ----------------------------------------------------------------
# Write the data out
# ----------------------------------------------------------------

#' Save the ALSPAC phenotype data onto a remote server
#'
#' @param outdat data.frame of ALSPAC phenotype data
#' @param filename name of file to store ALSPAC data (should be a .txt or .tsv file)
#' @param outpath path to output the data locally
#' @param remote_host hostname of remote server - needs to include username - i.e. USERNAME@REMOTE_SERVER
#' @param remote_user username for remote server
#' @param remote_dir directory on remote server to deposit the data
#' @param remote_pwd password to log onto the remote server
save_alsp_data <- function(outdat, filename, outpath, remote_host, remote_user, remote_dir, remote_pwd)
{
	ori_wd <- getwd()
	setwd(outpath)
	write.table(outdat, file = filename, 
				quote = F, col.names = T, row.names = F, sep = "\t")
	sftp_con <- sftp_connect(server = remote_host, 
				 username = remote_user, 
				 password = remote_pwd, 
				 folder = remote_dir, 
				 timeout = 500)
	sftp_upload(filename, fromfolder = ".", sftp_connection = sftp_con)
	message("Removing ALSPAC file from local directory")
	system(paste("rm", filename))
	setwd(ori_wd)
}

## FILL REST OF FIELDS IN MANUALLY AND DON'T SAVE THE ANSWERS
save_alsp_data(outdat = new_res, 
			   filename = outfile, 
			   outpath = outdir, 
			   remote_user = "", 
			   remote_host = "", 
			   remote_dir = "", 
			   remote_pwd = "")
