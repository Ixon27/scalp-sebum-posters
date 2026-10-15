.libPaths(c('~/R_libs', .libPaths()))
suppressPackageStartupMessages({
  library(ape)
  library(reshape2)
})

DATA <- 'data'
META <- 'data/metadata_sebum.txt'

# Read rarefied table
ft_rare <- read.table(file.path(DATA, 'feature_table_rare.tsv'),
                      header = TRUE, sep = '\t', skip = 1, comment.char = '',
                      check.names = FALSE)
colnames(ft_rare)[1] <- 'FeatureID'
cat('rarefied table dim:', dim(ft_rare), '\n')

ft_raw <- read.table(file.path(DATA, 'feature_table_raw.tsv'),
                     header = TRUE, sep = '\t', skip = 1, comment.char = '',
                     check.names = FALSE)
colnames(ft_raw)[1] <- 'FeatureID'
cat('raw table dim:', dim(ft_raw), '\n')

# Taxonomy
tax <- read.table(file.path(DATA, 'taxonomy/taxonomy.tsv'),
                  header = TRUE, sep = '\t', quote = '', comment.char = '')
colnames(tax)[1] <- 'FeatureID'

split_tax <- function(s) {
  parts <- strsplit(s, ';', fixed = TRUE)[[1]]
  parts <- trimws(parts)
  out <- setNames(rep(NA_character_, 7),
                  c('Domain','Phylum','Class','Order','Family','Genus','Species'))
  for (p in parts) {
    if (startsWith(p, 'd__')) out['Domain']  <- substring(p, 4)
    if (startsWith(p, 'p__')) out['Phylum']  <- substring(p, 4)
    if (startsWith(p, 'c__')) out['Class']   <- substring(p, 4)
    if (startsWith(p, 'o__')) out['Order']   <- substring(p, 4)
    if (startsWith(p, 'f__')) out['Family']  <- substring(p, 4)
    if (startsWith(p, 'g__')) out['Genus']   <- substring(p, 4)
    if (startsWith(p, 's__')) out['Species'] <- substring(p, 4)
  }
  out
}
tax_split <- as.data.frame(do.call(rbind, lapply(tax$Taxon, split_tax)),
                           stringsAsFactors = FALSE)
tax_split <- cbind(FeatureID = tax$FeatureID, tax_split)

# Fill blanks
for (col in c('Phylum','Class','Order','Family','Genus')) {
  tax_split[[col]][is.na(tax_split[[col]]) | tax_split[[col]] == ''] <- 'Unassigned'
}

# Aggregate to genus level (rarefied)
merge_genus <- function(ft, tax_split) {
  m <- merge(ft, tax_split[, c('FeatureID','Genus','Phylum')],
             by = 'FeatureID')
  samples <- setdiff(colnames(ft), 'FeatureID')
  agg <- aggregate(m[, samples], by = list(Genus = m$Genus), FUN = sum)
  # also keep dominant phylum per genus
  phy_map <- aggregate(m$Phylum, by = list(Genus = m$Genus),
                       FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
  list(table = agg, phy = setNames(phy_map$x, phy_map$Genus))
}

g_rare <- merge_genus(ft_rare, tax_split)
g_raw  <- merge_genus(ft_raw,  tax_split)
cat('genus-level rarefied dim:', dim(g_rare$table), '\n')

# Metadata
meta <- read.table(META, header = TRUE, sep = '\t', comment.char = '')
meta <- meta[meta$SampleID != '#q2:types', ]
for (col in c('Age','Burnish','Moisture','Sebum','Desquamation','Temperature','Redness')) {
  meta[[col]] <- as.numeric(meta[[col]])
}
rownames(meta) <- meta$SampleID
cat('metadata rows:', nrow(meta), '\n')

# Save objects
saveRDS(list(genus_rare = g_rare$table,
             genus_raw  = g_raw$table,
             phylum_map = g_rare$phy,
             meta       = meta,
             tax_split  = tax_split),
        file.path(DATA, 'prep.rds'))

# Relative abundance (rarefied)
samples <- setdiff(colnames(g_rare$table), 'Genus')
ra <- g_rare$table
ra[, samples] <- sweep(ra[, samples], 2, colSums(ra[, samples]), '/') * 100
write.table(ra, file.path(DATA, 'genus_rel_abund.tsv'), sep = '\t',
            quote = FALSE, row.names = FALSE)
cat('top 10 mean abundance:\n')
top10 <- ra[order(-rowMeans(ra[, samples])), ][1:10, ]
print(data.frame(Genus = top10$Genus, MeanRA = rowMeans(top10[, samples])))
cat('\nDone.\n')
