source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(patchwork)
})

DATA <- 'data'
OUT  <- 'output'

prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare
phy  <- prep$phylum_map

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
meta <- meta[samples, ]
ra <- gen
ra[, samples] <- sweep(ra[, samples], 2, colSums(ra[, samples]), '/') * 100
ra <- ra[ra$Genus != 'Unassigned', ]

top <- ra[order(-rowMeans(ra[, samples]))[1:25], ]
rownames(top) <- top$Genus
mat <- as.matrix(top[, samples])

# Row-wise log10(RA+0.01) then Z-score (across samples)
log_mat <- log10(mat + 0.01)
z_mat   <- t(scale(t(log_mat)))

# Row clustering
row_hc <- hclust(dist(z_mat), method = 'ward.D2')
row_order <- rownames(z_mat)[row_hc$order]

# Column ordering: NS then SS, within each by sebum
ns_s <- samples[meta[samples,'Group'] == 'NS']
ss_s <- samples[meta[samples,'Group'] == 'SS']
ns_s <- ns_s[order(meta[ns_s, 'Sebum'])]
ss_s <- ss_s[order(meta[ss_s, 'Sebum'])]
col_order <- c(ns_s, ss_s)

z_df <- reshape2::melt(z_mat, varnames = c('Genus','Sample'),
                       value.name = 'z')
z_df$Genus  <- factor(z_df$Genus,  levels = row_order)
z_df$Sample <- factor(z_df$Sample, levels = col_order)

# Phylum annotation
genus_phy <- phy[as.character(row_order)]
genus_phy[is.na(genus_phy) | genus_phy == ''] <- 'Other'
pal_phy <- c(Actinobacteriota = '#E79453', Firmicutes = '#00B050',
             Proteobacteria   = '#0070C0', Bacteroidota = '#E4C868',
             Cyanobacteria    = '#92D050', Deinococcota = '#BF9000',
             Chloroflexi      = '#C84E40', Other = 'grey70')
genus_phy[!(genus_phy %in% names(pal_phy))] <- 'Other'
phy_df <- data.frame(Genus = factor(row_order, levels = row_order),
                     Phylum = factor(genus_phy, levels = unique(genus_phy)),
                     x = 1)

# Sample annotations
samp_ann <- data.frame(
  Sample = factor(col_order, levels = col_order),
  Group  = meta[col_order, 'Group'],
  Sebum  = meta[col_order, 'Sebum']
)

# === Main heatmap ===
p_heat <- ggplot(z_df, aes(Sample, Genus, fill = z)) +
  geom_tile(color = 'white', linewidth = 0.15) +
  scale_fill_gradientn(colors = heatmap_div, na.value = 'grey90',
                       name = 'Z-score',
                       limits = c(-2, 2), oob = scales::squish) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = 'Top 25 genera abundance heatmap (log10 RA, Z-scored)',
       subtitle = 'Rows: Ward.D2 cluster on Z-scores  |  Cols: NS then SS, ordered by sebum') +
  theme_bw_pub(11) +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        axis.text.y  = element_text(face = 'italic', size = 9),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = 'gray30', linewidth = 0.4))

# Phylum side strip (right)
p_phy <- ggplot(phy_df, aes(x = x, y = Genus, fill = Phylum)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_manual(values = pal_phy, name = 'Phylum') +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'right',
        legend.title = element_text(face = 'bold'))

# Sample Group bar (top)
p_grp <- ggplot(samp_ann, aes(x = Sample, y = 1, fill = Group)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_manual(values = my_colors, name = 'Group') +
  scale_x_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'top',
        legend.title = element_text(face = 'bold'))

# Sample Sebum bar (top)
p_sebum <- ggplot(samp_ann, aes(x = Sample, y = 1, fill = Sebum)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_gradient(low = '#FFF2CC', high = '#C0000D', name = 'Sebum') +
  scale_x_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'top',
        legend.title = element_text(face = 'bold'))

# Row dendrogram (left) — base impl, no ggdendro dependency
hcdata <- {
  ord_data <- function(hc) {
    nL <- length(hc$labels)
    m  <- hc$merge
    leaf_x <- setNames(seq_len(nL), hc$labels[hc$order])
    segs <- list()
    node_y <- numeric(nrow(m))
    node_x <- numeric(nrow(m))
    for (i in seq_len(nrow(m))) {
      l <- m[i, 1]; r <- m[i, 2]
      xL <- if (l < 0) leaf_x[hc$labels[-l]] else node_x[l]
      xR <- if (r < 0) leaf_x[hc$labels[-r]] else node_x[r]
      yL <- if (l < 0) 0 else node_y[l]
      yR <- if (r < 0) 0 else node_y[r]
      h <- hc$height[i]
      segs[[length(segs) + 1]] <- data.frame(x = xL, xend = xL, y = yL, yend = h)
      segs[[length(segs) + 1]] <- data.frame(x = xR, xend = xR, y = yR, yend = h)
      segs[[length(segs) + 1]] <- data.frame(x = xL, xend = xR, y = h,  yend = h)
      node_x[i] <- (xL + xR) / 2
      node_y[i] <- h
    }
    do.call(rbind, segs)
  }
  ord_data(row_hc)
}

p_dend <- ggplot(hcdata, aes(x = y, y = x, xend = yend, yend = xend)) +
  geom_segment(linewidth = 0.4) +
  scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.5, length(row_order) + 0.5)) +
  theme_void()

# === Combine via patchwork ===
top_row    <- plot_spacer() + p_grp   + plot_spacer() + plot_layout(widths = c(0.12, 1, 0.10))
sebum_row  <- plot_spacer() + p_sebum + plot_spacer() + plot_layout(widths = c(0.12, 1, 0.10))
heat_row   <- p_dend       + p_heat   + p_phy         + plot_layout(widths = c(0.12, 1, 0.10))

final <- (top_row / sebum_row / heat_row) +
  plot_layout(heights = c(0.04, 0.04, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

save_fig(final, 'F2_anno_heatmap', OUT, 13, 8)
cat('Saved F2.\n')
