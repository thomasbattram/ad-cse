# ---------------------------------------------------------------
# Summarise the results of the EWAS
# ---------------------------------------------------------------

## Aim: To take results from the Joint location-and-scale score test EWAS and summarise them in plots + tables

## pkgs
library(tidyverse) # tidy code, data, plots
library(ewaff) # for plotting QQ plots
library(cowplot) # plotting on a grid
library(usefunc) # own package of useful functions
library(RColorBrewer)

## args
args <- commandArgs(trailingOnly = TRUE)
res_files <- args[1]
cse_summ_file <- args[2]
plotfile <- args[3]
summ_outfile <- args[4]

message("Arguments: ", args)

## args tests
# res_files <- paste0("results/ewas/jlssc-res-", c("cc", "no-cc", "no-covs"), ".tsv")
# cse_summ_file <- "results/summary-of-results.RData"
# plotfile <- "results/jlssc-ewas-man-qq.png"
# summ_outfile <- "results/jlssc-ewas-summary-of-results.RData"

res_files <- unlist(str_split(res_files, " "))

## data
mods <- gsub("jlssc-res-|.tsv", "", basename(res_files))
names(mods) <- res_files
all_res <- lapply(res_files, read_tsv)
names(all_res) <- mods
cse_summ <- new_load(cse_summ_file)

## cases and controls
sample_size <- list(N = getmode(all_res[[1]]$N), 
					N_case = getmode(all_res[[1]]$N_cases), 
					N_controls = getmode(all_res[[1]]$N_controls))

# ---------------------------------------------------------------
# Functions for setup of data and qq plots etc.
# ---------------------------------------------------------------

#' Make a QQ plot
#' 
#' @param res results from the EWAS
#' 
#' @return QQ plot made using the ewaff package
make_qq <- function(res)
{
	lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    lamb2 <- paste("lambda == ", comma(lamb))
	ewaff_qq <- ewaff.qq.plot(res$P, lambda.method = "none") + 
		theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lamb2, hjust = 0, vjust = 1, parse = TRUE) + 
        # labs(title = p_title) +
        theme(plot.title = element_blank()) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
		theme(text = element_text(size = 8), legend.position = "none")
}

#' Make a Manhattan plot
#' 
#' @param res results from the EWAS - this should be the combined results using the comb_res function
#' @param cpg_annotations genomic annotations for the CpG sites so the genomic position can be plotted
#' 
#' @return Manhattan plot
make_man <- function(res, cpg_annotations, sigp=1e-7, sugp=1e-5, highl = FALSE, cpgs_to_highl = "")
{
    res$name <- res$probeID
    res <- res %>%
        left_join(cpg_annotations)
    # to highlight
    if (highl) {
    	if (cpgs_to_highl == "") {
    		cpg_h <- res[res$P < sigp, ]$name
    	} else {
    		cpg_h <- cpgs_to_highl
    	}
    } else {
    	cpg_h <- ""
    }
    gg_man <- gg.manhattan(df = res, 
                           hlight = cpg_h, 
                           title = NULL, 
                           SNP = "name", 
                           CHR = "chr", 
                           BP = "position", 
                           P = "P", 
                           sig = sigp, 
                           sugg = sugp, 
                           lab = FALSE, 
                           colour = TRUE)
    gg_man <- gg_man + 
        theme(plot.title = element_blank(), text = element_text(size = 10),
        	  axis.text.x = element_text(angle = 90, size = 6))
    return(gg_man)
}

#' Remove points from a scatter plot where density is really high
#' @param x x-coordinates vector
#' @param y y-coordinates vector
#' @param resolution number of partitions for the x and y-dimensions.
#' @param max.per.cell maximum number of points per x-y partition.
#' @return index into the points that omits points from x-y partitions
#' so that each has at most \code{max.per.cell} points.
scatter_thinning <- function(x,y,resolution=100,max.per.cell=100) {
    x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
    y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
    z.cell <- x.cell * resolution + y.cell
    frequency.table <- table(z.cell)
    frequency <- rep(0,max(z.cell, na.rm=T))
    frequency[as.integer(names(frequency.table))] <- frequency.table
    f.cell <- frequency[z.cell]
    
    big.cells <- length(which(frequency > max.per.cell))
    sort(c(which(f.cell <= max.per.cell),
           sample(which(f.cell > max.per.cell),
                  size=big.cells * max.per.cell, replace=F)),
         decreasing=F)
}

get_lamda <- function(p)
{
	median(qchisq(p, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
}

# ---------------------------------------------------------------
# Plot QQs and Manhattans
# ---------------------------------------------------------------

## Make QQs
# dat_text <- all_lambda %>%
# 	mutate(label = paste("lambda == ", comma(lambda))) %>%
# 	dplyr::select(-lambda)

qq_list <- lapply(all_res, make_qq)
# ggsave(qq_outfile, plot = qq_p)

## Make Manhattans
annotation <- meffil::meffil.get.features("epic")
annotation <- annotation %>% 
    mutate(chr = gsub("chr", "", chromosome)) %>%
    mutate(chr = gsub("X", "23", chr)) %>% 
    mutate(chr = as.numeric(gsub("Y", "24", chr)))


man_list <- lapply(all_res, make_man, cpg_annotations = annotation)
# ggsave(man_outfile, plot = man_out)

## Plot together

## Need to do in three parts to add a title to each part
models <- c("Cell count adjusted", "Cell count unadjusted", "Completely unadjusted")
names(models) <- c("cc", "no-cc", "no-covs")
comb_plots <- lapply(seq_along(models), function(x) {
    model <- models[x]
    mod <- names(models)[x]
    ## Put QQ plot and Manhattan side-by-side
    plot <- cowplot::plot_grid(qq_list[[mod]], man_list[[mod]], 
                               nrow = 1, ncol = 2, 
                               rel_widths = c(1, 2))
    ## Make title
    title <- ggdraw() + 
        draw_label(
            model,
            fontface = 'bold',
            x = 0,
            hjust = 0
        ) +
        theme(
            # add margin on the left of the drawing canvas,
            #  so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
        )  
    ## Put title and plot together
    out <- cowplot::plot_grid(title, plot, ncol = 1, rel_heights = c(0.1, 1))
    return(out)
})

summ_plots_out <- cowplot::plot_grid(plotlist = comb_plots, nrow = 3, ncol = 1)
ggsave(plotfile, plot = summ_plots_out)

# ---------------------------------------------------------------
# Summarise results
# ---------------------------------------------------------------

## What we want at the end:

## Results at P<1x10-5 in one table:
# Model, CpG, Beta, SE, P, WB EWAS hit, Cell EWAS hit
## samplesizes
## Number of CpGs tested

all_res_tab <- bind_rows(all_res, .id = "Model")
# all_res_tab[all_res_tab$probeID == "cg18478105", "P"] <- c(1e-6, 2e-6, 8e-6)
ewaff_hits <- cse_summ$ewaff_hits
cse_hits <- cse_summ$initial_hits

n_hits <- length(unique(cse_hits$CpG))

omicwas_p_thresh <- 0.05 / n_hits
omic_rep_cpgs <- cse_hits %>%
	dplyr::filter(Method == "omicWAS") %>%
	dplyr::filter(P < omicwas_p_thresh) %>%
	pull("CpG")

cse_hits_sig <- cse_hits %>%
	dplyr::filter(CpG %in% omic_rep_cpgs)

sig_res <- all_res_tab %>%
	dplyr::filter(P < 1e-5) %>%
	mutate(`WB EWAS hit` = ifelse(probeID %in% ewaff_hits$CpG, TRUE, FALSE),
		   `Cell EWAS hit` = "FALSE") %>%
	dplyr::select(Model, CpG = probeID, Beta = BETA, SE, P, `WB EWAS hit`, `Cell EWAS hit`)

for (i in 1:nrow(sig_res)) {
	cpg <- sig_res[i, "CpG", drop = TRUE]
	if (cpg %in% cse_hits_sig$CpG) {
		cells <- cse_hits_sig %>%
			dplyr::filter(CpG %in% cpg) %>%
			pull(`Cell type`) %>%
			unique()

		sig_res[i, "Cell EWAS hit"] <- paste0(cells, collapse = ", ")
	}
}

out_summ <- list(res = sig_res, samplesizes = sample_size, n_cpg = nrow(all_res[[1]]))
save(out_summ, file = summ_outfile)

