# ----------------------------------------
# Compare EWAS results for methods for estimating cell-specific effects
# ----------------------------------------

## Aim: Compare methods for estimating cell specific-effects

## Date: 2022-01-23

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # for EWAS functions
library(aries) # for extracting aries samplesheet
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
meta_file <- args[1]
ewas_resdir <- args[2]
traits_file <- args[3]
out_file <- args[4]

# setwd("SCRATCH_SPACE")
# meta_file <- "data/metadata.tsv"
# ewas_resdir <- "results/ewas-res"
# traits_file <- "data/traits.txt"
# ## traits_file <- "data/test-traits.txt"
# hits_out_file <- "results/extracted-hits.RData"
# summ_out_file <- "results/summary.RData"
# rep_out_file <- "results/replication-data.RData"
# hits_plot_file <- "results/summary-of-ewas-hits.png"

## data
meta_file <- read_tsv(meta_file)
traits <- readLines(traits_file)
methods <- c("celldmc", "tca", "tcareg", "toast", "ewaff")

## structure:
## 1. Extract results at p<1e-5 from each EWAS (and save these results)
## 2. Plot a summary of hits per method across EWAS
## 3. Assess replication of hits across each EWAS
## 4. Assess differences between cs methods and ewaff
	## a. including checking if meta-analysis will recapitulate ewaff
## 5. Assess hits per cell type (plot along with cell type proportions?)

### SORT ARIES DATA
meth_file <- "../sims/data/aries_fom.RData"
aries_meth <- new_load(meth_file)
cpg_nams_and_rows <- tibble(rown = 1:nrow(aries_meth), cpg = rownames(aries_meth))
rm(aries_meth)

## Checking AHRR CpG present in smoking EWAS
# cpg_nams_and_rows[cpg_nams_and_rows$cpg == "cg05575921", ]
# all_res[["r6010"]]$ewaff %>% arrange(p)

## It's there!! woop woop

# ----------------------------------------
# Extract results at p<1e-5
# ----------------------------------------

#' Extract specific CpGs or CpGs associated with the trait of interest at a specific P value
#' 
#' @param res trait and method specific results output from the "ewas.R" script
#' @param cpg character vector of CpG sites to extract
#' @param p_threshold numeric. P value threshold - any associations with a P less than this will be extracted
#' @return tibble with columns: "CpG", "beta", "se", "p", "cell_type"
sort_res <- function(res, cpg = NULL, p_threshold = 1e-5) 
{
	sig_res <- res$p %>%
		as.data.frame() %>%
		rownames_to_column(var = "CpG") %>%
		pivot_longer(!CpG, names_to = "cell_type", values_to = "p") %>%
		dplyr::filter(p < p_threshold)

	if (!is.null(cpg)) {
		sig_res <- sig_res %>%
			dplyr::filter(CpG %in% cpg)
	}

	beta_df <- res$beta %>%
		as.data.frame() %>%
		rownames_to_column(var = "CpG") %>%
		pivot_longer(!CpG, names_to = "cell_type", values_to = "beta")
	
	out_res <- sig_res %>%
		left_join(beta_df) 

	if ("se" %in% names(res)) {
		se_df <- res$se %>%
			as.data.frame() %>%
			rownames_to_column(var = "CpG") %>%
			pivot_longer(!CpG, names_to = "cell_type", values_to = "se")
		out_res <- out_res %>%
			left_join(se_df)
	}

	out_res <- out_res %>%
		dplyr::select(CpG, beta, one_of("se"), p, cell_type) %>%
		as_tibble
	return(out_res)
}

# traits <- traits[1:50]
start_time <- proc.time()
all_res <- lapply(traits, function(trait) {
	print(trait)
	out <- lapply(methods, function(method) {
		print(method)
		res_file <- paste0(method, "/", trait, ".RData")
		res_path <- file.path(ewas_resdir, res_file)
		if (!file.exists(res_path)) {
			warning(res_path, " does not exist.")
			return(NULL)
		}
		ewas_res <- new_load(res_path)
		out <- sort_res(ewas_res, p_threshold = 1e-5)
		print(out)
		return(out)
	})
	names(out) <- methods
	return(out)
})
names(all_res) <- traits
time_taken <- proc.time() - start_time
time_taken # 4536.853 elapsed (~1.25hours)

save(all_res, file = hits_out_file)

# ----------------------------------------
# Plot a summary of hits per method
# ----------------------------------------

## bar chart with cell type groups and a group of "total" hits

## Get hits at conventional threshold from results
# all_res <- new_load(hits_out_file)
plot_res <- lapply(all_res, bind_rows, .id = "method") %>%
	bind_rows(.id = "trait") %>%
	dplyr::filter(p < 1e-7)

## remove instances where method couldn't run
all_files <- expand.grid(trait = traits, method = methods) %>%
	mutate(file = paste0(method, "/", trait, ".RData"), 
		   full_path = file.path(ewas_resdir, file))

completed_analyses <- sapply(traits, function(tra) {
	trait_files <- all_files[all_files$trait == tra, ]
	files_exist <- file.exists(trait_files$full_path)
	ifelse(all(files_exist), return(TRUE), return(FALSE))
})

bad_traits <- names(completed_analyses)[!completed_analyses]

plot_res <- dplyr::filter(plot_res, !trait %in% bad_traits)

plot_res_summ <- plot_res %>%
	group_by(cell_type, method) %>%
	summarize(hits = n())	

total_m_hits <- plot_res_summ %>%
	group_by(method) %>%
	summarize(hits = sum(hits)) %>%
	mutate(cell_type = "Combined")

comb_plot_res <- bind_rows(plot_res_summ, total_m_hits) %>%
	dplyr::filter(cell_type != ".")

ewaff_res_hits <- comb_plot_res[comb_plot_res$method == "ewaff", "hits", drop = T]
plot_text <- paste0("Number of associations across all conventional EWAS = ", ewaff_res_hits)
cols <- get_cb_palette()[1:length(methods)]
ct <- unique(plot_res_summ$cell_type)[unique(plot_res_summ$cell_type) != "."]
x_axis_labs <- c(ct, "Combined")
summ_plot <- ggplot(comb_plot_res, aes(x = cell_type, y = hits, group = method, fill = method)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	scale_fill_manual(values = cols) + 
	scale_x_discrete(limits = x_axis_labs) + 
	annotate(geom = "text", x = 4, y = 7500, label = plot_text, colour = "red", size = 3) + 
	labs(x = "Cell type", y = "Associations at P < 1x10-7") +
	theme_bw()

ggsave(hits_plot_file, plot = summ_plot)

## Summary across traits (see if any driving things)
hit_summ <- plot_res %>%
	group_by(cell_type, method, trait) %>%
	summarize(n_hits = n()) %>%
	group_by(cell_type, method) %>%
	summarize(max = max(n_hits), min = min(n_hits), median = median(n_hits), mean = mean(n_hits), 
			  total = sum(n_hits))

# ----------------------------------------
# Replication across methods 
# ----------------------------------------

trait_meth_combos <- plot_res %>%
	dplyr::select(trait, method) %>% 
	distinct()
all_rep_res <- lapply(1:nrow(trait_meth_combos), function(i) {
	print(i)
	tmc <- trait_meth_combos[i, ]
	res_of_int <- plot_res %>%
		dplyr::filter(trait == tmc$trait, method == tmc$method)
	cpg_of_int <- res_of_int$CpG
	files_to_check <- all_files %>%
			dplyr::filter(trait == tmc$trait, method != tmc$method)
	other_res <- map_dfr(1:nrow(files_to_check), function(x) {
		f2c <- files_to_check[x, ]
		ewas_res <- new_load(f2c$full_path)
		if (f2c$method == "ewaff") {
			ewas_res <- lapply(ewas_res, function(er) { names(er) <- cpg_nams_and_rows$cpg; return(er) })
		}
		out <- sort_res(ewas_res, cpg = cpg_of_int, p_threshold = 1)
		out$method <- f2c$method
		return(out)
	})
	out <- other_res %>%
		dplyr::filter(cell_type %in% c(res_of_int$cell_type, ".")) %>%
		mutate(trait = tmc$trait) %>%
		bind_rows(res_of_int)
	return(out)
})

save(all_rep_res, file = rep_out_file)

## WRITE CODE TO ASSESS REPLICATION PROPERLY HERE!! 
# 1. direction of effect
# 2. p value replication at p<0.05/N_traits
# 3. heterogeneity of effects between CellDMC and TCA

# ----------------------------------------
# Meta-analysis of res
# ----------------------------------------

## Does meta-analysing res give same results as conventional EWAS?
# use CellDMC and TCA for this! 

## NEED TO ADD IN TCA RES AND FIGURE OUT HOW TO SORT DATA PROPERLY

traits_to_test <- 10
all_cpgs <- cpg_nams_and_rows$cpg
set.seed(2)
test_cpgs <- sample(all_cpgs, 100)
good_traits <- traits[!traits %in% bad_traits]
lapply(1:traits_to_test, function(x) {
	test_t <- good_traits[x]
	file_res <- all_files %>%
		dplyr::filter(trait == test_t, method %in% c("celldmc", "tca", "ewaff"))
	ewas_res <- new_load(file_res$full_path[1])
	ewaff_res <- new_load(file_res$full_path[3])
	ewaff_res <- lapply(ewaff_res, function(er) { names(er) <- cpg_nams_and_rows$cpg; return(er) })
	out <- sort_res(ewas_res, p_threshold = 1, cpg = test_cpgs)
	ewaff_out <- sort_res(ewaff_res, p_threshold = 1, cpg = test_cpgs)
	celldmc_test_res <- map_dfr(test_cpgs, function(cpg) {
		in_data <- out %>%
			dplyr::filter(CpG == cpg)
		m_gen <- metagen(TE = beta, 
					 seTE = se, 
					 studlab = cell_type, 
					 data = in_data, 
					 fixed = TRUE, 
					 random = FALSE)	
		tibble(CpG = cpg, beta = m_gen$TE.fixed)
	})
	temp <- left_join(celldmc_test_res, ewaff_out, by = c("CpG" = "CpG"))
	cor(temp$beta.x, temp$beta.y)
})

m.gen <- metagen(TE = TE,
                 seTE = seTE,
                 studlab = Author,
                 data = ThirdWave,
                 sm = "SMD",
                 fixed = TRUE,
                 random = FALSE)


## To do:
# 1. Get normal EWAS results for BMI from EWAS Catalog res
# 2. Make plots
# 	a. Normal EWAS "hits" vs. their results in specific cells
# 	b. Cell-specific "hits" that aren't in the normal EWAS
# 	c. Complete comparison of effect sizes? - maybe a heatmap with all methods + normal EWAS

## 2a.

get_hits <- function(res, p_threshold) dplyr::filter(res, p.value < p_threshold) %>% arrange(p.value)

names(all_res)

celldmc_res <- bind_rows(map(all_res, "celldmc"), .id = "trait")
tca_res <- bind_rows(map(all_res, "tca"), .id = "trait")
tcareg_res <- bind_rows(map(all_res, "tcareg"), .id = "trait")
toast_res <- bind_rows(map(all_res, "toast"), .id = "trait")
ewaff_res <- bind_rows(map(all_res, "ewaff"), .id = "trait")

lapply(methods, function(method) {
	temp <- bind_rows(map(all_res, method), .id = "trait")
	table(temp$cell_type)
})

## SHIT ON IT!! UNSURE WHAT CPGS ARE IN EWAFF RES!!!

x <- 1e-7/50

tca_res %>%
	dplyr::filter(CpG %in% celldmc_res$CpG)
celldmc_res %>%
	dplyr::filter(CpG %in% tca_res$CpG)


celldmc_res %>% dplyr::filter(p < x)
tca_res %>% dplyr::filter(p < x)
tcareg_res %>% dplyr::filter(p < x)
toast_res %>% dplyr::filter(p < x)

## BMI
ori_res <- get_hits(original_ewas[["BMI"]], 1e-7) %>%
	dplyr::select(CpG = probeID, beta = estimate, se, p = p.value) %>%
	mutate(cell_type = "bulk")
cs_res <- lapply(methods, function(method) {
		res_file <- paste0(method, "/", "fm1ms111.RData")
		res_path <- file.path(ewas_resdir, res_file)
		ewas_res <- new_load(res_path)
		out <- sort_res(ewas_res, cpg = ori_res$CpG, p_threshold = 1)
		return(out)
})
names(cs_res) <- methods
cs_res <- cs_res[names(cs_res) != "omicwas"]
cs_res$celldmc <- dplyr::select(cs_res$celldmc, -se)

cs_res_df <- bind_rows(cs_res, .id = "method")

plot_dat <- ori_res %>%
	dplyr::select(-se) %>%
	mutate(method = "bulk") %>%
	bind_rows(cs_res_df) %>%
	mutate(effect_dir = ifelse(sign(beta) == -1, "neg", "pos"), logp = -log10(p))

cb_pal <- get_cb_palette()
cols <- cb_pal[1:2]
names(cols) <- c("neg", "pos")
abcg1_bar <- ggplot(plot_dat[plot_dat$method != "bulk",], aes(x = method, y = logp, fill = as.factor(effect_dir))) + 
	geom_bar(stat="identity") +
	geom_hline(yintercept = plot_dat[plot_dat$method == "bulk", "logp", drop=T], colour = cols[2]) +  
	scale_fill_manual(name = "DoE", values = cols) + 
	labs(y = "-log10(P)") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle=90)) + 
	facet_grid(~ cell_type)

ggsave("abcg1-comp.png", abcg1_bar)

## Smoking
smok_ori_res <- get_hits(original_ewas[["Smoking"]], 1e-7) %>%
	dplyr::select(CpG = probeID, beta = estimate, se, p = p.value) %>%
	mutate(cell_type = "bulk")
smok_cs_res <- lapply(methods, function(method) {
		res_file <- paste0(method, "/", "r6010.RData")
		res_path <- file.path(ewas_resdir, res_file)
		ewas_res <- new_load(res_path)
		out <- sort_res(ewas_res, cpg = smok_ori_res$CpG, p_threshold = 1)
		return(out)
})
names(smok_cs_res) <- methods
smok_cs_res <- smok_cs_res[names(smok_cs_res) != "omicwas"]

ahrr_res <- lapply(smok_cs_res, function(x) {dplyr::filter(x, CpG == "cg05951221")})

ahrr_res$celldmc <- dplyr::select(ahrr_res$celldmc, -se)

ahrr_res_df <- bind_rows(ahrr_res, .id = "method")

plot_dat <- smok_ori_res %>%
	dplyr::filter(CpG == "cg05951221") %>%
	dplyr::select(-se) %>%
	mutate(method = "bulk") %>%
	bind_rows(ahrr_res_df) %>%
	mutate(effect_dir = ifelse(sign(beta) == -1, "neg", "pos"), logp = -log10(p))


ahhr_bar <- ggplot(plot_dat[plot_dat$method != "bulk",], aes(x = method, y = logp, fill = as.factor(effect_dir))) + 
	geom_bar(stat="identity") +
	geom_hline(yintercept = plot_dat[plot_dat$method == "bulk", "logp", drop=T], colour = cols[1]) +  
	scale_fill_manual(name = "DoE", values = cols) + 
	labs(y = "-log10(P)") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle=90)) + 
	facet_grid(~ cell_type)

ggsave("ahrr-comp.png", ahhr_bar)


### tabs of both
smok_res_list <- lapply(all_res$Smoking, function(x){ 
	if(nrow(x) == 0) return(NULL)
	return(x)
})

smok_ori_all_res <- original_ewas$Smoking %>%

smok_res_df <- smok_res_list %>%
	bind_rows(.id = "method") %>%
	left_join(original_ewas$Smoking, by = c("CpG" = "probeID")) %>%
	dplyr::select(-Details)

write.table(smok_res_df, file = "smok_res_df_out.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

cpg_of_int <- smok_res_df$

res_file <- paste0(method, "/", trait, ".RData")
res_path <- file.path(ewas_resdir, "tca/r6010.RData")
ewas_res <- new_load(res_path)

temp_out <- sort_res(ewas_res, cpg = cpg_of_int, p_threshold = 1)
write.table(temp_out, file = "smok-tca-cg23248424.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

all_res3 <- lapply(methods, function(method) {
	print(method)
	res_file <- paste0(method, "/r6010.RData")
	res_path <- file.path(ewas_resdir, res_file)
	ewas_res <- new_load(res_path)
	out <- sort_res(ewas_res, p_threshold = 1, cpg = cpg_of_int)
	return(out)
})
names(all_res3) <- methods
all_res3 <- all_res3[names(all_res3) != "omicwas"]
all_res3$celldmc <- dplyr::select(all_res3$celldmc, -se)


all_res3_df <- bind_rows(all_res3, .id = "method")
write.table(all_res3_df, file = "smok-allmeth-cg23248424.tsv", sep = "\t", quote = F, row.names = F, col.names = T)










