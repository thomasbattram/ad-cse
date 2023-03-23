# ---------------------------------------------------------------
# Compare extended blood cell counts with cell counts generated with previous reference
# ---------------------------------------------------------------

## Compare cell counts generated using the cell counts estimated using meffil and using EpiDISH.

## pkgs
library(tidyverse) # tidy code and data
library(cowplot) # plots

## args
args <- commandArgs(trailingOnly = T)
meffil_cc_file <- args[1] # cell counts created using meffil
epidish_cc_file <- args[2] # cell counts created using EpiDISH
boxplot_outfile <- args[3] # boxplot comparing the cell counts + methods
scatter_outfile <- args[4] # scatter plot comparing the methods

## manual args
meffil_cc_file <- "data/extended-blood-celltypes-idats.tsv" # cell counts created using meffil 
epidish_cc_file <- "data/extended-blood-celltypes-epidish.tsv" # cell counts created using EpiDISH
boxplot_outfile <- "cc-comparison-plot.png" # boxplot comparing the cell counts + methods
scatter_outfile <- "cc-comparison-stats.tsv" # scatter plot comparing the methods

## data
meffil_cc <- read_tsv(meffil_cc_file)
epi_cc <- read_tsv(epidish_cc_file)

# ---------------------------------------------------------------
# Make boxplot plot
# ---------------------------------------------------------------

common_samples <- intersect(meffil_cc$IID, epi_rpc$IID)

all_dat <- bind_rows(list(meffil = meffil_cc, epidish = epi_cc), .id = "Algorithm") %>%
	dplyr::filter(IID %in% common_samples)

## Format for plotting
plot_dat <- all_dat %>%
	pivot_longer(cols = -c(Algorithm, IID), names_to = "Cell type", values_to = "Proportion")

outplot <- ggplot(plot_dat, aes(x = `Cell type`, y = Proportion, fill = Algorithm)) + 
	theme_bw() + 
	geom_boxplot() +
	theme(legend.position = "bottom")

ggsave(boxplot_outfile, plot = outplot, width = 12)

# ---------------------------------------------------------------
# Scatter plots
# ---------------------------------------------------------------

celltypes <- unique(plot_dat[["Cell type"]])

scat_list <- lapply(celltypes, function(ct) {
	p_dat <- all_dat %>%
		dplyr::select(Algorithm, IID, one_of(ct)) %>%
		pivot_wider(names_from = Algorithm, values_from = one_of(ct))

	scatplot <- ggplot(p_dat, aes_string(x = meffil, y = epidish)) + 
		geom_point() + 
		geom_abline(intercept = 0, slope = 1, colour = "red") + 
		theme_bw()
	return(scatplot)
})

out_scatter <- cowplot::plot_grid(plotlist = scat_list, nrow = 3)
ggsave(scatter_outfile, out_scatter)
