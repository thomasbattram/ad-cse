# ----------------------------------------------------------------
# Extract ARIES phenotype data for comparison of EWAS
# ----------------------------------------------------------------

## Aim: To extract eczema data from ALSPAC participants, using two different definitions of AD to compare

## This is for data at age 7

## pkgs
library(alspac) # extract data
library(tidyverse) # tidy code and data
library(openxlsx) # write out data as an excel spreadsheet
library(sftp)
library(labelled) # label data
library(usefunc) # personal package of useful functions

## args
# aries_samplesfile <- "data/aries-samples-f7.tsv" # THIS NEEDS TO BE DOWNLOADED FROM BC4/BP1 AND DELETED AFTER ANALYSES ARE COMPLETE
# samplesizes_file <- "data/ad-def-samplesizes.xlsx"
# outdir <- "data"
# outfile <- "ad-data-f7.tsv"

aries_samples <- read_tsv(aries_samplesfile)

# ----------------------------------------------------------------
# sort alspac
# ----------------------------------------------------------------

alspac_data_dir <- "/Volumes/Data/"
setDataDir(alspac_data_dir)

data(current)
# data(useful)

## RUN THIS TO UPDATE THE DICTIONARIES
# current <- createDictionary("Current", name="current")
# useful <- createDictionary("Useful_data", name="useful")

# ----------------------------------------------------------------
# Find variables
# ----------------------------------------------------------------

## Here we are defining Childhood AD - any reported eczema that has seen a doctor (or not)
## Variables are taken from the PACE analyses of Childhood AD
ecz_vars <- c("kq035", "kr042", "ks1042", "kv1060", "kv1070") 
current %>%
	dplyr::filter(name %in% ecz_vars) %>%
	dplyr::select(name, lab)

## The date at which these variables were ascertained can be found using the ALSPAC data dictionary
## 		- Timepoints should be near when DNAm was taken

## Timepoints for variables above:
# kq035, A3w: CH Had Eczema In Past Year - timepoint = 6y 9m
# kr042, A3w: Child had eczema in past 12 months - timepoint = 7y 7m
# ks1042, A3w: Child had eczema in past year - timepoint = 8y 7m
# kv1060, A4u: Child had eczema in past 12 months - timepoint = 10y 8m
# kv1070, Doctor stated that child has asthma or eczema - timepoint = 10y 8m 

ecz_timepoints <- tibble(`ALSPAC variable` = ecz_vars, 
						 timepoint = c("6y 9m", "7y 7m", "8y 7m", "10y 8m", "10y 8m"))

# ----------------------------------------------------------------
# Extract variables and labels
# ----------------------------------------------------------------

new_current <- current %>%
	dplyr::filter(name %in% ecz_vars)

## extraction of data
result <- extractVars(new_current)

## extraction of labels
val_labs <- usefunc::extract_alspac_labels(new_current, alsp_dir = alspac_data_dir)
unique_labs <- unique(unlist(val_labs))

# ----------------------------------------------------------------
# Define AD
# ----------------------------------------------------------------

aries_res <- aries_samples %>%
	left_join(result) %>%
	dplyr::select(aln, alnqlet, qlet, Sample_Name, all_of(ecz_vars)) %>%
	mutate(ad_dr = case_when(kq035 == 1  |
							  kr042 == 1  | 
							  ks1042 == 1 | 
							  kv1060 == 1 |
							  kv1070 == 2 | 
							  kv1070 == 3 ~ 1,
							  kq035 == 3  |
							  kr042 == 3  | 
							  ks1042 == 3 | 
							  kv1060 == 3 |
							  kv1070 == 4 ~ 0), 
		   ad_both = case_when(kq035 == 1  |
		   					  kq035 == 2  |
							  kr042 == 1  | 
							  kr042 == 2  | 
							  ks1042 == 1 | 
							  ks1042 == 2 | 
							  kv1060 == 1 |
							  kv1060 == 2 |
							  kv1070 == 2 | 
							  kv1070 == 3 ~ 1,
							  kq035 == 3  |
							  kr042 == 3  | 
							  ks1042 == 3 | 
							  kv1060 == 3 |
							  kv1070 == 4 ~ 0))

# write out sample sizes and definitions
## TAB:
# | alspac variable | label | notes | timepoint | N | cases | controls |

cases <- list(ad_dr = sum(aries_res$ad_dr == 1, na.rm=T), 
			  ad_both = sum(aries_res$ad_both == 1, na.rm=T)) 

controls <- list(ad_dr = sum(aries_res$ad_dr == 0, na.rm=T), 
			  	 ad_both = sum(aries_res$ad_both == 0, na.rm=T))

sample_summary <- tibble(`AD variable` = c("ad_doc", "ad_all"), 
						 notes = c("Cases = doctor diagnosed only - reported eczema with no doctor = set to missing",
						   		   "Cases = Any report of AD - doctor diagnosed or not."), 
						 N = c(cases$ad_dr + controls$ad_dr, cases$ad_both + controls$ad_both),
						 cases = c(cases$ad_dr, cases$ad_both), 
						 controls = c(controls$ad_dr, controls$ad_both))

## Add to other samplesizes file
wb <- loadWorkbook(samplesizes_file)
addWorksheet(wb, "f7") 
writeData(wb, "f7", sample_summary)
saveWorkbook(wb, samplesizes_file, overwrite = TRUE)

## Using "ad_both" for now
clean_res <- aries_res %>%
	mutate(ad = ad_both) %>%
	dplyr::select(aln, alnqlet, qlet, ad) %>%
	na.omit() %>%
	as_tibble()

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
				 folder = remote_dir)
	sftp_upload(filename, fromfolder = ".", sftp_connection = sftp_con)
	message("Removing ALSPAC file from local directory")
	system(paste("rm", filename))
	setwd(ori_wd)
}

## FILL REST OF FIELDS IN MANUALLY AND DON'T SAVE THE ANSWERS
save_alsp_data(outdat = clean_res, 
			   filename = outfile, 
			   outpath = outdir, 
			   remote_user = "", 
			   remote_host = "", 
			   remote_dir = "", 
			   remote_pwd = "")

system(paste("rm", aries_samplesfile))
