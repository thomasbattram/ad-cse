# ----------------------------------------------------------------
# Extract ARIES phenotype data for comparison of EWAS
# ----------------------------------------------------------------

## Aim: To extract eczema data from ALSPAC participants, using two different definitions of AD to compare

## pkgs
library(alspac) # extract data
library(tidyverse) # tidy code and data
library(openxlsx) # write out data as an excel spreadsheet
library(sftp)
library(labelled) # label data
library(usefunc) # personal package of useful functions

## args
# aries_samplesfile <- "data/aries-samples-15up.tsv" # THIS NEEDS TO BE DOWNLOADED FROM BC4/BP1 AND DELETED AFTER ANALYSES ARE COMPLETE
# samplesizes_file <- "data/ad-def-samplesizes.xlsx"
# outdir <- "data"
# outfile <- "ad-data.tsv"

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

ecz_vars <- findVars(c("eczema", "atopic dermatitis"))
ecz_vars %>%
	dplyr::filter(grepl("Child", cat3))

## The date at which these variables were ascertained can be found using the ALSPAC data dictionary
## 		- Timepoints should be near when DNAm was taken
##		- Could also do a "teenage eczema" phenotype

## Based on this, here are some potential variables of interest:
# ccs5023, I3d: YP has had eczema in the past 12 months - timepoint = 16y 6m
# cct4055, C19: Respondent ever had eczema - timepoint = 18y
# tc6110, F11: Study teenager has had eczema - timepoint = 16y 6m
# ta1030, A2u:  Teenager had eczema in past 12 months - timepoint = 13y 1m
# tb1060, A4u: Child had eczema in past 12 months - timepoint = 13y 10m 
# tb1070, A5: Doctor ever diagnosed asthma/eczema - timepoint = 13y 10m
# tb1122, A9c: Child ever had eczema - timepoint = 13y 10m
# tb1170, A10e1: Number of days child taken off school for eczema/itchy rash in past 12 mo - timepoint = 13y 10m
# tb1171, A10e2: Guess at number of days child taken off school for eczema/itchy rash in p - timepoint = 13y 10m
# tb2210, B1v: Medicine/pills/drops/ointment used by child for eczema in past 12 months - timepoint = 13y 10m
# tb2213, B1v: Frequency child used medicine/pills/drops/ointment for eczema in past 12 mo - timepoint = 13y 10m

## variables to check the sample size of 
ecz_vars_to_check <- c("ccs5023", "cct4055", "tc6110", "ta1030", "tb1060", "tb1070", "tb2210")

ecz_timepoints <- tibble(`ALSPAC variable` = ecz_vars_to_check, 
						 timepoint = c("16y 6m", "18y", "16y 6m", "13y 1m", "13y 10m", "13y 10m", "13y 10m"))

# ----------------------------------------------------------------
# Extract variables and labels
# ----------------------------------------------------------------

new_current <- current %>%
	dplyr::filter(name %in% ecz_vars_to_check)

## extraction of data
result <- extractVars(new_current)

## extraction of labels
val_labs <- usefunc::extract_alspac_labels(new_current, alsp_dir = alspac_data_dir)
unique_labs <- unique(unlist(val_labs))

# ----------------------------------------------------------------
# Check samplesizes for each variable
# ----------------------------------------------------------------

aries_res <- aries_samples %>%
	left_join(result) %>%
	dplyr::select(aln, alnqlet, qlet, all_of(ecz_vars_to_check))

## Turn values into missing (looks like negatives should be set to missing)
aries_res[aries_res < 0] <- NA

## Replace values with case/control
## NOTE: for cases where there is an option for "Yes, but didn't see a doctor"
# 		 they will be set to missing and not controls when just looking at doctor diagnosis only
aries_res <- aries_res %>%
	mutate(ta1030_dr = case_when(ta1030 == 1 ~ "case", 
								 ta1030 == 3 ~ "control"), 
		   ta1030_both = case_when(ta1030 == 1 | ta1030 == 2 ~ "case", 
								   ta1030 == 3 ~ "control"),
		   tb1060_dr = case_when(tb1060 == 1 ~ "case", 
		   						 tb1060 == 3 ~ "control"), 
		   tb1060_both = case_when(tb1060 == 1 | tb1060 == 2 ~ "case", 
								   tb1060 == 3 ~ "control"),
		   tb1070 = case_when(tb1070 == 2 | tb1070 == 3 ~ "case", 
		   					  tb1070 == 1 | tb1070 == 4 ~ "control"), 
		   tb2210 = case_when(tb2210 == 1 ~ "case", 
		   					  tb2210 == 0 ~ "control"), 
		   tc6110 = case_when(tc6110 == 1 ~ "case", 
		   					  tc6110 == 2 ~ "control"), 
		   ccs5023_dr = case_when(ccs5023 == 1 ~ "case", 
		   						  ccs5023 == 3 ~ "control"), 
		   ccs5023_both = case_when(ccs5023 == 1 | ccs5023 == 2 ~ "case", 
		   						  	ccs5023 == 3 ~ "control"),
		   cct4055 = case_when(cct4055 == 1 ~ "case", 
		   					   cct4055 == 2 ~ "control")
			) 

aries_res <- aries_res %>%
	mutate(ad_teen_dr = case_when(ta1030_dr == "case" | tb1060_dr == "case" |
								  tb1070 == "case" | tb2210 == "case" | 
								  ccs5023_dr == "case" ~ "case", 
								  ta1030_dr == "control" | tb1060_dr == "control" |
								  tb1070 == "control" | tb2210 == "control" | 
								  ccs5023_dr == "control" ~ "control"), 
		   ad_late_teen_dr = case_when(ccs5023_dr == "case" ~ "case", 
								  	   ccs5023_dr == "control" ~ "control"), 
		   ad_teen_both = case_when(ta1030_both == "case" | tb1060_both == "case" |
								  	tb1070 == "case" | tb2210 == "case" | 
								  	ccs5023_both == "case" ~ "case", 
								  	ta1030_both == "control" | tb1060_both == "control" |
								  	tb1070 == "control" | tb2210 == "control" | 
								  	ccs5023_both == "control" ~ "control"),
		   ad_late_teen_both = case_when(ccs5023_both == "case" ~ "case", 
								  		 ccs5023_both == "control" ~ "control"), 								  	 
		   ad_ever = case_when(tc6110 == "case" | cct4055 == "case" ~ "case", 
		   					   tc6110 == "control" | cct4055 == "control" ~ "control")
		   )

## NOTE: ad_ever is from mid-late teen questionnaires - could go earlier and up sample size

## ADD IN A COMBINED PHENOTYPE!!

lapply(aries_res[, -c(1:3)], function(x) {table(x)})

# write out sample sizes and definitions
## TAB:
# | alspac variable | label | notes | timepoint | N | cases | controls |

vars <- c("cct4055", "tc6110", "ccs5023_dr", "ccs5023_both", "ta1030_dr",
		  "ta1030_both", "tb1060_dr", "tb1060_both", "tb1070", "tb2210", 
		  "ad_teen_dr", "ad_late_teen_dr", "ad_teen_both", "ad_late_teen_both", 
		  "ad_ever")
sample_summary <- map_dfr(vars, function(var) {
	if (grepl("ad_", var)) {
		lab <- ""
		tp <- ""
		cur_var <- var
	} else {
		cur_var <- gsub("_.*", "", var)
		filt_cur <- new_current[new_current$name == cur_var, ]	
		lab <- filt_cur$lab
		tp <- ecz_timepoints[ecz_timepoints[["ALSPAC variable"]] == cur_var, "timepoint", drop=T]
	}
	case_n <- sum(aries_res[[var]] == "case", na.rm=T)
	control_n <- sum(aries_res[[var]] == "control", na.rm=T)
	if (grepl("dr", var)) notes <- "Cases are those that went to a doctor"
	if (grepl("both", var)) notes <- "Cases include those that did and did not go to a doctor"
	if (!grepl("both|dr", var)) notes <- ""
	tibble(`ALSPAC variable` = cur_var, 
		   label = lab, 
		   notes = notes, 
		   timepoint = tp,
		   N = case_n + control_n, 
		   cases = case_n, 
		   controls = control_n)
})

write.xlsx(sample_summary, file = samplesizes_file)

## Using "ad_teen_both" for now
clean_res <- aries_res %>%
	mutate(ad = case_when(ad_teen_both == "case" ~ 1, 
						  ad_teen_both == "control" ~ 0)) %>%
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
