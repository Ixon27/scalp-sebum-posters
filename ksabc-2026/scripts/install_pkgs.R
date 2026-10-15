.libPaths(c('~/R_libs', .libPaths()))
options(repos = c(CRAN = 'https://cran.rstudio.com/'), Ncpus = max(1, parallel::detectCores()-2))
cran_pkgs <- c('dplyr','tidyr','ggpubr','patchwork','corrplot','ggridges','tibble','readr','stringr','scales','BiocManager')
to_inst <- setdiff(cran_pkgs, rownames(installed.packages()))
if(length(to_inst)) install.packages(to_inst, lib='~/R_libs')

library(BiocManager, lib.loc='~/R_libs')
bioc_pkgs <- c('ggtree','ggtreeExtra','treeio','ComplexHeatmap','microbiome')
to_inst_b <- setdiff(bioc_pkgs, rownames(installed.packages()))
if(length(to_inst_b)) BiocManager::install(to_inst_b, ask=FALSE, update=FALSE, lib='~/R_libs')

cat('\n=== final check ===\n')
all_p <- c(cran_pkgs, bioc_pkgs)
for(p in all_p) cat(sprintf('%-18s : %s\n', p, requireNamespace(p, quietly=TRUE)))
