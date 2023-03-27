pkgs <- list(cran=c("devtools", "tidyverse", "BiocManager", "matrixStats", "cowplot", 
                    "cluster", "readxl", "SmartSVA", "TCA"),
             bioc=c("EpiDISH", "sva"),
             git=c("https://github.com/perishky/meffil", 
                   "https://github.com/perishky/ewaff", 
                   "https://github.com/thomasbattram/usefunc", 
                   "https://github.com/jrs95/jlst"))

for (pkg in pkgs$cran) {
  cat("R package:", pkg, "\n")
  installed <- installed.packages()[,"Package"]
  if (!pkg %in% installed)
     install.packages(pkg)
}

for (pkg in pkgs$bioc) {
  cat("R package:", pkg, "\n")
  installed <- installed.packages()[,"Package"]
  if (!pkg %in% installed)
    BiocManager::install(pkg)
}

for (url in pkgs$git) {
  installed <- installed.packages()[,"Package"]
  pkg <- basename(url)
  cat("R package:", pkg, "\n")
  if (!pkg %in% installed)
    devtools::install_github(url)
}    

