# ---------------------------------------------------------------
# Compare extended blood cell counts with cell counts generated with previous reference
# ---------------------------------------------------------------

## Compare cell counts generated using the "blood gse167998" and "blood gse35069 complete" references.

## pkgs
library(tidyverse) # tidy code and data
library(aries) # easy access to the aries data on bc4/bp1
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = T)
aries_dir <- args[1]
new_cc_file <- args[2]
vio_plot_file <- args[3]
cor_stats_file <- args[4]

## manual args
aries_dir <- ""
new_cc_file <- "../data-extraction-and-qc/data/extended-blood-celltypes-epic-15up-idats.tsv"
vio_plot_file <- "cc-comparison-plot.png"
cor_stats_file <- "cc-comparison-stats.tsv"

## data
new_cc <- read_tsv(new_cc_file)
aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
old_cc <- aries$cell.counts[["blood-gse35069-complete"]] 

# ---------------------------------------------------------------
# Put celltypes into categories to enable comparisons
# ---------------------------------------------------------------

## celltype broad categories
colnames(new_cc)
colnames(old_cc)

celltype_cats <- tibble(celltype = c("Treg", "CD4nv", "CD4mem", "CD4T", 
									 "CD8nv", "CD8mem", "CD8T", 
									 "Bnv", "Bmem", "Bcell", 
									 "Bas", "Neu", "Eos",
									 "Mono", "NK"),
						broad_celltype = c(rep("CD4+ T", 4), 
										   rep("CD8+ T", 3), 
										   rep("B cell", 3), 
										   rep("Granulocyte", 3), 
										   rep("Other", 2))
						)


# ---------------------------------------------------------------
# Make violin/boxplot plot
# ---------------------------------------------------------------

new_cc_plot <- new_cc %>%
	pivot_longer(!Sample_Name, names_to = "celltype", values_to = "value") %>%
	mutate(reference = "gse167998")	%>%
	left_join(celltype_cats) %>%
	dplyr::filter(Sample_Name %in% rownames(old_cc))

old_cc_plot <- as.data.frame(old_cc)
old_cc_plot$Sample_Name <- rownames(old_cc)

old_cc_plot <- old_cc_plot %>%
	pivot_longer(!Sample_Name, names_to = "celltype", values_to = "value") %>%
	mutate(reference = "gse35069")	%>%
	left_join(celltype_cats)

plot_res <- bind_rows(new_cc_plot, old_cc_plot)


## plot_res should look like this:
## broad_celltype | celltype | reference | value
## 	T cell        | T reg    | gse167998 | ...
## ...


outplot <- ggplot(plot_res, aes(x = celltype, y = value, fill = reference)) + 
	# geom_violin() + 
	theme_bw() + 
	geom_boxplot() +
	# geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) + 
	facet_grid(. ~broad_celltype, scales = "free") + 
	theme(legend.position = "bottom")

ggsave(vio_plot_file, plot = outplot, width = 12)

# ---------------------------------------------------------------
# Comparison stats
# ---------------------------------------------------------------

## Make sure tables are comparable
new_cc_stats <- new_cc %>%
	dplyr::filter(Sample_Name %in% rownames(old_cc)) %>%
	arrange(Sample_Name)

old_cc_stats <- as.data.frame(old_cc) %>%
	rownames_to_column(var = "Sample_Name") %>%
	dplyr::filter(Sample_Name %in% new_cc_stats$Sample_Name) %>%
	arrange(Sample_Name)

stopifnot(all(old_cc_stats$Sample_Name == new_cc_stats$Sample_Name))

## Make comparison across all cell types derived from old reference 
celltypes <- c("Mono", "NK", "Bcell", "CD4T", "CD8T")
ct_comp_stats_list <- lapply(celltypes, function(ct) {
	broad_ct <- celltype_cats %>%
		dplyr::filter(celltype == ct) %>%
		pull(broad_celltype)
	old_res <- old_cc_stats[[ct]]
	if (broad_ct == "other") {
		new_res <- new_cc_stats[[ct]]
	} else {
		ccc <- celltype_cats %>%
			dplyr::filter(broad_celltype == broad_ct)
		new_cts <- ccc$celltype[ccc$celltype != ct]
		new_res <- rowSums(new_cc_stats[, which(colnames(new_cc_stats) %in% new_cts)])
	}
	p_out_cor <- cor(new_res, old_res, method = "pearson")
	s_out_cor <- cor(new_res, old_res, method = "spearman")
	t_out_diff <- t.test(new_res, old_res, paired = TRUE, alternative = "two.sided")
	w_out_diff <- wilcox.test(new_res, old_res, paired = TRUE, alternative = "two.sided")	
	new_cc_sumstats <- extract_sum_stats(new_res)
	out_stats <- extract_sum_stats(old_cc_stats[[ct]]) %>%
		bind_rows(new_cc_sumstats) %>%
		mutate(celltype = ct, reference = c("gse35069", "gse167998"), 
			   pearson_r = p_out_cor, spearman_r = s_out_cor, 
			   ttest_p = t_out_diff$p.value, wilcoxtest_p = w_out_diff$p.value) %>%
		dplyr::select(celltype, reference, everything())
	return(out_stats)
})
ct_comp_stats <- bind_rows(ct_comp_stats_list)

write.table(ct_comp_stats, file = cor_stats_file, col.names = T, row.names = F, sep = "\t", quote = F)

## Correlation of identical cell types

## Correlation of combined cell types

## summary of both 



## Mono and NK cells should be most consistent
## Bmem + Bnv should equal B cells
## Bas + Eos + Neu (new) should equal Neu + Eos (old)
## Unsure how T regs change things for T cells...


### ALL CELLTYPES

## T cells
## T Reg +
## CD4mem + 
## CD4nv + 
## CD4 -
## CD8mem + 
## CD8nv +
## CD8 -

## B cells
## Bmem + 
## Bnv + 
## B -

## Granulocytes
## Neu +-
## Eos +-
## Bas +

## Other
## Mono +-
## NK +-