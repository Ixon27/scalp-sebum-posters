source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(Hmisc)
  library(reshape2)
})

DATA <- 'data'
OUT  <- 'output'

prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare
phy  <- prep$phylum_map

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
ra <- gen
ra[, samples] <- sweep(ra[, samples], 2, colSums(ra[, samples]), '/') * 100

# Drop 'Unassigned' for plot clarity
ra <- ra[ra$Genus != 'Unassigned', ]

# Top 30 by mean abundance
ord <- order(-rowMeans(ra[, samples]))
top30 <- ra[ord[1:30], ]
rownames(top30) <- top30$Genus
top30_mat <- t(as.matrix(top30[, samples]))  # samples x genera

clin_vars <- c('Age','Burnish','Moisture','Sebum','Desquamation','Temperature','Redness')
clin_mat <- as.matrix(meta[samples, clin_vars])

# Spearman correlation (genera as columns, clinical as rows)
mat_all <- cbind(top30_mat, clin_mat)
res <- rcorr(mat_all, type = 'spearman')
r_mat <- res$r[clin_vars, colnames(top30_mat)]
p_mat <- res$P[clin_vars, colnames(top30_mat)]

# To long
df <- melt(r_mat, varnames = c('Clinical','Genus'), value.name = 'r')
df$p <- melt(p_mat)$value
df$r[is.na(df$r)] <- 0
df$p[is.na(df$p)] <- 1
df$star <- ''
df$star[df$p < 0.05]  <- '*'
df$star[df$p < 0.01]  <- '**'
df$star[df$p < 0.001] <- '***'

# Order genera by overall |r| sum (or keep abundance order)
df$Genus    <- factor(df$Genus, levels = colnames(top30_mat))
df$Clinical <- factor(df$Clinical, levels = rev(clin_vars))

# Phylum strip for genera (top annotation)
genus_phy <- phy[colnames(top30_mat)]
genus_phy[is.na(genus_phy) | genus_phy == ''] <- 'Other'
phy_df <- data.frame(Genus = factor(colnames(top30_mat), levels = colnames(top30_mat)),
                     Phylum = genus_phy)

pal_phy <- c(Actinobacteriota = '#E76F51', Firmicutes = '#2A9D8F',
             Proteobacteriota = '#264653', Proteobacteria = '#264653',
             Bacteroidota = '#F4A261', Cyanobacteria = '#8AB17D',
             Deinococcota = '#9B5DE5', Chloroflexi = '#BDB2FF',
             Other = '#999999', Unassigned = '#cccccc')
phy_df$Phylum_col <- ifelse(phy_df$Phylum %in% names(pal_phy),
                            as.character(phy_df$Phylum), 'Other')

p_main <- ggplot(df, aes(x = Genus, y = Clinical, fill = r)) +
  geom_tile(color = 'white', linewidth = 0.3) +
  geom_text(aes(label = star), size = 4, vjust = 0.78, color = 'black') +
  scale_fill_gradientn(colors = heatmap_div, limits = c(-1, 1),
                       na.value = '#eeeeee', name = 'Spearman r',
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = NULL, y = NULL,
       title = 'Genus-clinical correlations (top 30 genera)',
       subtitle = '* p<0.05  ** p<0.01  *** p<0.001') +
  theme_bw_pub(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        panel.grid = element_blank())

p_strip <- ggplot(phy_df, aes(x = Genus, y = 1, fill = Phylum_col)) +
  geom_tile() +
  scale_fill_manual(values = pal_phy, name = 'Phylum') +
  theme_void() +
  theme(legend.position = 'right',
        plot.margin = margin(0, 0, 0, 0))

save_fig(p_main, 'F4_spearman_heatmap', OUT, 11, 4.2)

cat('Saved F4 to', OUT, '\n')

# Quick sanity: top correlations involving Sebum
sebum_row <- df[df$Clinical == 'Sebum', ]
sebum_row <- sebum_row[order(-abs(sebum_row$r)), ]
cat('Top |r| with Sebum:\n')
print(head(sebum_row[, c('Genus','r','p','star')], 10))
