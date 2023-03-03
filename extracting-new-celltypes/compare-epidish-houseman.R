# -----------------------------------------------------------------
# Comapre epidish and houseman cell count estimates
# -----------------------------------------------------------------

# srun --job-name "InteractiveJob" --account=smed001801 --partition=veryshort --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=2:00:00 --mem=50GB --pty bash

## pkgs
library(EpiDISH)
library(aries)

new_cc_file <- "../data-extraction-and-qc/data/extended-blood-celltypes-epic-15up-idats.tsv"
aries_dir <- "/user/work/ms13525/aries"

## calc cell counts using epidish
cc <- readr::read_tsv(new_cc_file)
aries <- aries.select(aries_dir, time.point = "15up", featureset = "epic")
meth <- aries.methylation(aries)
meth <- meth[, colnames(meth) %in% cc$Sample_Name]
rm(cc)
data(cent12CT.m)
out.l <- epidish(beta.m = meth, ref.m = cent12CT.m, method = "RPC")
head(out.l$estF)

rm(meth)

library(tidyverse)
cc <- read_tsv(new_cc_file)
epi_cc <- as_tibble(out.l$estF, rownames = "Sample_Name")

## data needs to be in format:

## Model | Cell type | value


colnames(cc)
## combine data
common_samples <- intersect(cc$Sample_Name, epi_cc$Sample_Name)
all_dat <- epi_cc %>%
	rename(CD4nv = CD4Tnv, Bas = Baso, CD4mem = CD4Tmem, CD8mem = CD8Tmem, 
		   CD8nv = CD8Tnv) %>%
	bind_rows(cc, .id = "Algorithm") %>%
	mutate(Algorithm = ifelse(Algorithm == 1, "EpiDISH", "Houseman")) %>%
	dplyr::filter(Sample_Name %in% common_samples)



## Format for plotting
plot_dat <- all_dat %>%
	pivot_longer(cols = -c(Algorithm, Sample_Name), names_to = "Cell type", values_to = "Proportion")

outplot <- ggplot(plot_dat, aes(x = `Cell type`, y = Proportion, fill = Algorithm)) + 
	# geom_violin() + 
	theme_bw() + 
	geom_boxplot() +
	# geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) + 
	theme(legend.position = "bottom")

plot_file <- "epidish-houseman-boxplot.png"
ggsave(plot_file, plot = outplot, width = 12)
