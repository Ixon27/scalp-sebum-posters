source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggnewscale)
  library(ggrepel)
  library(patchwork)
})

# Drop-in for Hmisc::rcorr (Spearman): identical r (Pearson on ranks) and p
# (t-approximation: t = r*sqrt((n-2)/(1-r^2))) -> network edges unchanged.
rcorr <- function(m, type = 'spearman') {
  R <- cor(m, method = 'spearman')
  n <- nrow(m)
  t <- R * sqrt((n - 2) / (1 - R^2))
  P <- 2 * pt(-abs(t), n - 2)
  diag(P) <- NA
  list(r = R, P = P)
}

DATA <- 'data'
OUT  <- 'output'
prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare
phy  <- prep$phylum_map

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
ra <- gen
ra[, samples] <- sweep(ra[, samples], 2, colSums(ra[, samples]), '/') * 100
ra <- ra[ra$Genus != 'Unassigned', ]

ord <- order(-rowMeans(ra[, samples]))
top <- ra[ord[1:30], ]
rownames(top) <- top$Genus

ns_samp <- samples[meta[samples,'Group'] == 'NS']
ss_samp <- samples[meta[samples,'Group'] == 'SS']
mat_ns <- t(as.matrix(top[, ns_samp]))
mat_ss <- t(as.matrix(top[, ss_samp]))

R_CUT <- 0.6
P_CUT <- 0.05

corr_edges <- function(m) {
  cc <- rcorr(m, type = 'spearman')
  r <- cc$r; p <- cc$P
  el <- which(upper.tri(r), arr.ind = TRUE)
  df <- data.frame(from = rownames(r)[el[,1]], to = rownames(r)[el[,2]],
                   r = r[upper.tri(r)], p = p[upper.tri(p)])
  df <- df[!is.na(df$r) & !is.na(df$p), ]
  df[abs(df$r) >= R_CUT & df$p < P_CUT, ]
}

e_ns <- corr_edges(mat_ns); cat('NS edges:', nrow(e_ns), '\n')
e_ss <- corr_edges(mat_ss); cat('SS edges:', nrow(e_ss), '\n')

# Shared layout from union graph (comparable panels)
all_edges <- unique(rbind(e_ns[, c('from','to')], e_ss[, c('from','to')]))
all_nodes <- top$Genus
g_union <- graph_from_data_frame(all_edges, vertices = data.frame(name = all_nodes), directed = FALSE)
set.seed(42)
layout_xy <- layout_with_fr(g_union, niter = 2000)
rownames(layout_xy) <- V(g_union)$name

make_g <- function(edges, nodes)
  graph_from_data_frame(edges[, c('from','to','r')],
                        vertices = data.frame(name = nodes), directed = FALSE)
g_ns <- make_g(e_ns, all_nodes); g_ss <- make_g(e_ss, all_nodes)
deg_ns <- degree(g_ns); deg_ss <- degree(g_ss)

# Muted, poster-cohesive phylum palette (avoids navy/red, reserved for NS/SS)
pal_phy <- c(Actinobacteriota = '#C7913B', Firmicutes = '#2E7E8C',
             Proteobacteria   = '#5B7DB1', Bacteroidota = '#8A8A46',
             Cyanobacteria    = '#6FA287', Deinococcota = '#9C6B4F',
             Chloroflexi      = '#7E6E9C', Other = 'grey72')

make_edge_df <- function(e) {
  d <- e[, c('from','to','r')]
  d$x    <- layout_xy[d$from, 1]; d$y    <- layout_xy[d$from, 2]
  d$xend <- layout_xy[d$to,   1]; d$yend <- layout_xy[d$to,   2]
  d$sign <- ifelse(d$r > 0, 'Positive', 'Negative')
  d
}
make_node_df <- function(deg) data.frame(
  Genus = all_nodes, x = layout_xy[all_nodes,1], y = layout_xy[all_nodes,2],
  Phylum = ifelse(all_nodes %in% names(phy) & phy[all_nodes] %in% names(pal_phy),
                  phy[all_nodes], 'Other'),
  MeanRA = rowMeans(top[, samples]), deg = deg[all_nodes])

edge_sign_pal <- c(Positive = 'grey60', Negative = '#08519c')
edge_lty_pal  <- c(Positive = 'solid',  Negative = '22')

build_panel <- function(edf, ndf, ttl, sub, title_col = 'black') {
  lab <- ndf %>% filter(deg >= 4) %>% arrange(desc(deg)) %>% head(8)
  ggplot() +
    geom_segment(data = edf,
                 aes(x = x, y = y, xend = xend, yend = yend,
                     color = sign, linetype = sign, linewidth = abs(r)),
                 alpha = 0.65, lineend = 'round') +
    scale_color_manual(values = edge_sign_pal, name = 'Co-occurrence') +
    scale_linetype_manual(values = edge_lty_pal, name = 'Co-occurrence') +
    scale_linewidth_continuous(range = c(0.3, 1.4), guide = 'none') +
    new_scale_color() +
    # white halo so each node reads clearly above the edge web
    geom_point(data = ndf, aes(x = x, y = y, size = MeanRA),
               shape = 21, fill = 'white', color = 'white', stroke = 2.6) +
    geom_point(data = ndf, aes(x = x, y = y, size = MeanRA, fill = Phylum),
               shape = 21, color = 'gray25', stroke = 0.5) +
    scale_fill_manual(values = pal_phy, name = 'Phylum') +
    scale_size_continuous(range = c(2.5, 14), name = 'Mean RA (%)', breaks = c(1, 5, 20, 50)) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    geom_text_repel(data = lab, aes(x = x, y = y, label = Genus),
                    fontface = 'italic', size = 3.7, color = 'black', seed = 7,
                    box.padding = 0.6, point.padding = 0.4,
                    force = 5, force_pull = 0.2, max.overlaps = Inf,
                    min.segment.length = 0,
                    segment.color = 'gray60', segment.size = 0.3) +
    coord_fixed(clip = 'off') +
    labs(title = ttl, subtitle = sub) +
    theme_void(base_size = 13) +
    theme(text = element_text(family = 'sans'),
          plot.title = element_text(face = 'bold', size = 18, hjust = 0.5, color = title_col),
          plot.subtitle = element_text(size = 12.5, color = 'black', hjust = 0.5),
          plot.margin = margin(6, 6, 6, 6),
          plot.background = element_rect(fill = 'white', color = NA),
          legend.position = 'right',
          legend.title = element_text(face = 'bold', size = 13.5),
          legend.text = element_text(size = 12))
}

# panel metric strings
pmet <- function(g, e) sprintf('%d edges  ·  mean degree %.1f  ·  %d isolated genera',
                               nrow(e), 2 * nrow(e) / vcount(g), sum(degree(g) == 0))

edf_ns <- make_edge_df(e_ns); ndf_ns <- make_node_df(deg_ns)
edf_ss <- make_edge_df(e_ss); ndf_ss <- make_node_df(deg_ss)

p_ns <- build_panel(edf_ns, ndf_ns, 'NS  —  normal scalp', pmet(g_ns, e_ns),
                    title_col = '#0040B8')
p_ss <- build_panel(edf_ss, ndf_ss, 'SS  —  seborrheic scalp', pmet(g_ss, e_ss),
                    title_col = '#C0000D')

p_net <- (p_ns | p_ss) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'Co-occurrence network collapse in seborrheic scalp',
    theme = theme(plot.title = element_text(face = 'bold', size = 22),
                  plot.subtitle = element_text(size = 12.5, color = 'black'),
                  plot.background = element_rect(fill = 'white', color = NA))) &
  theme(legend.position = 'right')

save_fig(p_net, 'F3a_diff_network_themed', OUT, 14, 7)

cat('Saved F3a themed (muted phylum palette + navy/red panel titles).\n')
