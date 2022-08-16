# ---------------------------------------------------------------
# Summarise the results of the EWAS
# ---------------------------------------------------------------

## Aim: To take results from the cell-specific effects EWAS and summarise them in plots + tables

## pkgs
library(tidyverse) # tidy code, data, plots
library(ewaff) # for plotting QQ plots
library(cowplot) # plotting on a grid
library(usefunc) # own package of useful functions
library(RColorBrewer)

## args
args <- commandArgs(trailingOnly = TRUE)
celldmc_file <- args[1]
tca_file <- args[2]
omicwas_file <- args[3]
celldmc_tca_hits_file <- args[4]
ewaff_file <- args[5]
man_outfile <- args[6]
qq_outfile <- args[7]
summ_outfile <- args[8]

## args tests
# celldmc_file <- "results/ewas/celldmc-res.RData"
# tca_file <- "results/ewas/tca-res.RData"
# omicwas_file <- "results/ewas/omicwas-res.RData"
# celldmc_tca_hits_file <- "results/celldmc-tca-hits.RData"
# ewaff_file <- "results/ewas/ewaff-res.tsv"
# man_outfile <- "results/celldmc-tca-manhattans.png"
# qq_outfile <- "results/celldmc-tca-qqs.png"
# summ_outfile <- "results/summary-of-results.RData"

## data
celldmc <- new_load(celldmc_file)
tca <- new_load(tca_file)
omicwas <- new_load(omicwas_file)
hits <- new_load(celldmc_tca_hits_file)
hits_note <- hits$note
hits <- hits$res
ewaff_res <- read_tsv(ewaff_file)

# ---------------------------------------------------------------
# Functions for setup of data and qq plots etc.
# ---------------------------------------------------------------

#' Combine the summary statistics from each cell type
#' 
#' @param res results from the EWAS
#' @param celltype the cell type to combine results for
#' 
#' @return tibble of results with the summary statistics and "CpG" as columns 
comb_res <- function(res, celltype) {
	new_res <- map(res, as_tibble, rownames = "CpG")
	out <- map_dfc(new_res, celltype) %>%
		mutate(CpG = new_res$p$CpG) %>%
		dplyr::select(CpG, everything())
	return(out)
}


#' Make a QQ plot
#' 
#' @param res results from the EWAS - this should be the combined results using the comb_res function
#' 
#' @return QQ plot made using the ewaff package
make_qq <- function(res, p_title)
{
	lamb <- median(qchisq(res$p, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    lamb2 <- paste("lambda == ", comma(lamb))
	ewaff_qq <- ewaff.qq.plot(res$p, lambda.method = "none") + 
		theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lamb2, hjust = 0, vjust = 1, parse = TRUE) + 
        labs(title = p_title) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
		theme(text = element_text(size = 8), legend.position = "none")
}

#' Make a Manhattan plot
#' 
#' @param res results from the EWAS - this should be the combined results using the comb_res function
#' @param cpg_annotations genomic annotations for the CpG sites so the genomic position can be plotted
#' 
#' @return Manhattan plot
make_man <- function(res, cpg_annotations)
{
    res$name <- res$CpG
    res <- res %>%
        left_join(cpg_annotations)
    # to highlight
    cpg_h <- res[res$p < 1e-7, ]$name
    gg_man <- gg.manhattan(df = res, 
                           hlight = cpg_h, 
                           title = NULL, 
                           SNP = "name", 
                           CHR = "chr", 
                           BP = "position", 
                           P = "p", 
                           sig = 1e-7, 
                           sugg = 1e-5, 
                           lab = TRUE, 
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


setup_for_qq <- function(res, celltype)
{
	stats <- res[[celltype]] %>%
		mutate(expected = -log(sort(ppoints(p), decreasing = T), 10), 
		   	   observed = -log(sort(p, decreasing=T), 10)) %>%
		dplyr::select(expected, observed)
	selection.idx <- scatter_thinning(stats$observed, stats$expected,
        							  resolution = 100, max.per.cell = 100)
	out <- stats[selection.idx, ]
	return(out)
}

get_lamda <- function(p)
{
	median(qchisq(p, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
}

# ---------------------------------------------------------------
# Plot QQs and Manhattans
# ---------------------------------------------------------------

celltypes <- names(celldmc$beta)
celldmc_res <- lapply(celltypes, comb_res, res = celldmc)
names(celldmc_res) <- celltypes

tca_res <- lapply(celltypes, comb_res, res = tca)
names(tca_res) <- celltypes

## Plot QQs
tca_qq_res <- sapply(celltypes, setup_for_qq, res = tca_res, simplify = FALSE, USE.NAMES = TRUE) %>%
	bind_rows(.id = "celltype") %>%
	mutate(method = "tca") %>%
	dplyr::select(method, celltype, expected, observed)

all_qq_res <- sapply(celltypes, setup_for_qq, res = celldmc_res, simplify = FALSE, USE.NAMES = TRUE) %>%
	bind_rows(.id = "celltype") %>%
	mutate(method = "celldmc") %>%
	dplyr::select(method, celltype, expected, observed) %>%
	bind_rows(tca_qq_res)

all_lambda <- map_dfr(celltypes, function(x) {
	tca_l <- tibble(method = "tca", celltype = x, lambda = get_lamda(tca_res[[x]]$p))
	celldmc_l <- tibble(method = "celldmc", celltype = x, lambda = get_lamda(celldmc_res[[x]]$p))
	bind_rows(tca_l, celldmc_l)
})

dat_text <- all_lambda %>%
	mutate(label = paste("lambda == ", comma(lambda))) %>%
	dplyr::select(-lambda)

qq_p <- ggplot(all_qq_res, aes(x = expected, y = observed)) + 
	geom_abline(intercept = 0, slope = 1, colour = "black") + 
	geom_point() + 
	theme_bw() + 
    facet_grid(celltype ~ method) + 
    geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1, parse = TRUE)

ggsave(qq_outfile, plot = qq_p)

## Plot Manhattans
annotation <- meffil::meffil.get.features("epic")
annotation <- annotation %>% 
    mutate(chr = gsub("chr", "", chromosome)) %>%
    mutate(chr = gsub("X", "23", chr)) %>% 
    mutate(chr = as.numeric(gsub("Y", "24", chr)))

tca_man_res <- tca_res %>%
	bind_rows(.id = "celltype") %>%
	mutate(method = "tca") %>%
	dplyr::select(method, celltype, CpG, p)

all_man_res <- celldmc_res %>%
	bind_rows(.id = "celltype") %>%
	mutate(method = "celldmc") %>%
	dplyr::select(method, celltype, CpG, p) %>%
	bind_rows(tca_man_res)

man_out <- make_man(all_man_res, annotation) + 
	facet_grid(celltype ~ method)

ggsave(man_outfile, plot = man_out)

# WHAT IS A CELL TYPE?

# Epigenomic and transcriptomic analyses define core cell types, genes and targetable mechanisms for kidney disease

# ---------------------------------------------------------------
# Summarise results
# ---------------------------------------------------------------

## Here are the results that replicated across all 3 methods and the cell 

## FOR EACH RESULT:

## TAB1: Results for just hits in their cell types across each method
## | Cell type | Method | CpG | Beta | SE | P |

omicwas_res <- lapply(celltypes, comb_res, res = omicwas)
names(omicwas_res) <- celltypes

all_hit_res <- lapply(1:length(hits), function(x) {
	res <- hits[[x]]
	res_nam <- names(hits)[x]
	res <- res[res$replicated, ]
	if (nrow(res) == 0) return(NULL)
	celltype <- unlist(str_split(res_nam, "_"))[3]
	tca_filt_res <- tca_res[[celltype]] %>%
		dplyr::filter(CpG %in% res$CpG)
	celldmc_filt_res <- celldmc_res[[celltype]] %>%
		dplyr::filter(CpG %in% res$CpG) %>%
		dplyr::select(-se)
	omicwas_filt_res <- omicwas_res[[celltype]] %>%
		dplyr::filter(CpG %in% res$CpG)
	omicwas_rep <- omicwas_filt_res[omicwas_filt_res$p < 0.05, "CpG", drop=T]
	out <- bind_rows(list(TCA = tca_filt_res, CellDMC = celldmc_filt_res, omicWAS = omicwas_filt_res), 
					 .id = "Method") %>%
		   mutate(`Cell type` = celltype) %>%
		   dplyr::select(Method, `Cell type`, CpG, Beta = beta, P = p)
		   # dplyr::filter(CpG %in% omicwas_rep)
	return(out)
})
all_hit_res <- bind_rows(all_hit_res) %>%
	arrange(`Cell type`, CpG) %>%
	distinct()

sig_hits <- lapply(celltypes, function(ct) {
	all_hit_res %>%
		dplyr::filter(`Cell type` == ct) %>%
		dplyr::filter(Method == "omicWAS", P < 0.05) %>%
		pull(CpG)
})
names(sig_hits) <- celltypes
sig_hits <- unlist(sig_hits)

## FOR TESTING
# sig_hits <- "cg05173528"

ewaff_res <- ewaff_res %>%
	dplyr::select(CpG = probeID, Beta = BETA, SE = SE, P)

methods <- c("celldmc", "tca", "omicwas")
all_out <- lapply(methods, function(method) {
	res <- get(paste0(method, "_res"))
	map_dfr(sig_hits, function(sh) {
		map_dfr(celltypes, function(ct) {
			res[[ct]] %>%
				dplyr::filter(CpG == sh) %>%
				mutate(`Cell type` = ct) %>%
				dplyr::select(`Cell type`, CpG, Beta = beta, P = p)
		})
	})
})
names(all_out) <- methods

all_out$ewaff <- ewaff_res[ewaff_res$CpG %in% sig_hits, ]

## TO DO!!! - maybe think about best way when actually have some results...

## TAB2 - TABX: Results across other cell types for each CpG 
## METHOD ONE
## | Cell type | CpG | Beta | SE | P |
## | ...
## METHOD TWO
## | Cell type | CpG | Beta | SE | P |
## | ...
## METHOD THREE
## | Cell type | CpG | Beta | SE | P |
## | ...

summ_out <- list(initial_hits = all_hit_res, all_res = all_out)
save(summ_out, file = summ_outfile)








