.libPaths(c('~/R_libs', .libPaths()))
options(repos = c(CRAN = 'https://cran.rstudio.com/'), Ncpus = max(1, parallel::detectCores()-2))
library(BiocManager)
to_inst <- setdiff(c('treeio','ggtree','ggtreeExtra'), rownames(installed.packages()))
if(length(to_inst)) BiocManager::install(to_inst, ask=FALSE, update=FALSE,
                                          lib='~/R_libs',
                                          INSTALL_opts='--no-build-vignettes')
for(p in c('treeio','ggtree','ggtreeExtra'))
  cat(sprintf('%-15s : %s\n', p, requireNamespace(p, quietly=TRUE)))
