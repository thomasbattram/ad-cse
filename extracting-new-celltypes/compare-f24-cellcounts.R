# ---------------------------------------------------------------
# Compare extended blood cell counts with cell counts generated with previous reference
# ---------------------------------------------------------------

## Compare cell counts generated using the "blood gse167998" and "blood gse35069 complete" references.

## pkgs
library(aries)
library(tidyverse) # tidy code and data
library(cowplot) # arranging plots on grid nicely
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = T)
aries_dir <- args[1]
dir_cc_file <- args[2]
new_cc_file <- args[3]
vio_plot_file <- args[4]
cor_stats_file <- args[5]
hist_plot_file <- args[6]

## manual args
aries_dir <- ""
dir_cc_file <- "f24-direct-cellcounts.tsv"
new_cc_file <- "../../data-extraction-and-qc/data/extended-blood-celltypes-epic-f24-idats.tsv"
vio_plot_file <- "f24-cc-comparison-plot.png"
cor_stats_file <- "f24-cc-comparison-stats.tsv"
hist_plot_file <- "f24-direct-cellcounts-hist.png"

## data
dir_cc <- read_tsv(dir_cc_file)
new_cc <- read_tsv(new_cc_file)

# ---------------------------------------------------------------
# plot histogram of directly measured cell counts
# ---------------------------------------------------------------

celltypes <- c("Neutrophils_F24", "Lymphocytes_F24", "Monocytes_F24", "Eosinophils_F24",
		  	   "Basophils_F24")

hist_list <- lapply(celltypes, function(ct) {
	ggplot(dir_cc, aes_string(x = ct)) + 
		geom_histogram(fill = "blue", colour = "black") + 
		labs(x = paste(gsub("_F24", "", ct), "(10^9/L)")) +
		theme_bw()
})

histplot <- plot_grid(plotlist = hist_list, nrow = 3)
ggsave(hist_plot_file, plot = histplot)

# ---------------------------------------------------------------
# calculate proportion of cell types
# ---------------------------------------------------------------

## Limit directly measured cell counts to those in the F24 ARIES data
aries <- aries.select(aries_dir, time.point = "F24")
samples <- aries$samples %>%
	dplyr::select(alnqlet, Sample_Name)

dir_cc_filt <- dir_cc %>%
	na.omit() %>%
	left_join(samples) %>%
	dplyr::filter(Sample_Name %in% new_cc$Sample_Name)

## Calculate total number of cells per L
dir_cc_filt$total <- rowSums(dir_cc_filt[, c(celltypes)])

## Match estimated cellcount names to directly measured cell counts
colnames(new_cc)
names(celltypes) <- c("Neu", "Lymphocytes", "Mono", "Eos", "Bas")

## Calculate proportions
dir_cc_temp <- tibble(Sample_Name = dir_cc_filt$Sample_Name)
dir_cc_prop_list <- lapply(1:length(celltypes), function(x) {
	ct <- celltypes[x]
	prop <- dir_cc_filt[[ct]] / dir_cc_filt$total
	out <- dir_cc_temp
	out[[names(celltypes)[x]]] <- prop
	return(out)
})

dir_cc_prop <- reduce(dir_cc_prop_list, left_join)
rowSums(dir_cc_prop[, names(celltypes)])

## Get plot data ready!
new_cc_plot <- new_cc %>%
	pivot_longer(!Sample_Name, names_to = "celltype", values_to = "value") %>%
	mutate(method = "estimated")	%>%
	dplyr::filter(Sample_Name %in% dir_cc_prop$Sample_Name)

dir_cc_plot <- dir_cc_prop %>%
	pivot_longer(!Sample_Name, names_to = "celltype", values_to = "value") %>%
	mutate(method = "directly measured")	%>%
	dplyr::filter(Sample_Name %in% new_cc_plot$Sample_Name)

cts <- intersect(unique(new_cc_plot$celltype), unique(dir_cc_plot$celltype))

plot_res <- bind_rows(new_cc_plot, dir_cc_plot) %>%
	dplyr::filter(celltype %in% cts)

outplot <- ggplot(plot_res, aes(x = celltype, y = value, fill = method)) + 
	# geom_violin() + 
	theme_bw() + 
	geom_boxplot() +
	# geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) + 
	theme(legend.position = "bottom")

ggsave(vio_plot_file, plot = outplot, width = 12)

new_cc_filt <- new_cc[new_cc$Sample_Name %in% dir_cc_prop$Sample_Name,] %>%
	arrange(Sample_Name)
dir_cc_prop <- dir_cc_prop %>%
	arrange(Sample_Name)
lapply(cts, function(ct) {
	cor(new_cc_filt[[ct]], dir_cc_prop[[ct]])
})
