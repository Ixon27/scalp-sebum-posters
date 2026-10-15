source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(igraph)
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
OUT  <- 'output/v2'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Rebuild NS / SS networks exactly as in 03_F3 (top 30 genera)
# ============================================================
prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
ra <- gen
ra[, samples] <- sweep(ra[, samples], 2, colSums(ra[, samples]), '/') * 100
ra <- ra[ra$Genus != 'Unassigned', ]
ord <- order(-rowMeans(ra[, samples]))
top <- ra[ord[1:30], ]; rownames(top) <- top$Genus

ns_samp <- samples[meta[samples,'Group'] == 'NS']
ss_samp <- samples[meta[samples,'Group'] == 'SS']
mat_ns <- t(as.matrix(top[, ns_samp]))
mat_ss <- t(as.matrix(top[, ss_samp]))

R_CUT <- 0.6; P_CUT <- 0.05
corr_edges <- function(m) {
  cc <- rcorr(m, type = 'spearman'); r <- cc$r; p <- cc$P
  el <- which(upper.tri(r), arr.ind = TRUE)
  df <- data.frame(from = rownames(r)[el[,1]], to = rownames(r)[el[,2]],
                   r = r[upper.tri(r)], p = p[upper.tri(p)])
  df <- df[!is.na(df$r) & !is.na(df$p), ]
  df[abs(df$r) >= R_CUT & df$p < P_CUT, ]
}
e_ns <- corr_edges(mat_ns); e_ss <- corr_edges(mat_ss)
all_nodes <- top$Genus
make_g <- function(edges, nodes)
  graph_from_data_frame(edges[, c('from','to','r')],
                        vertices = data.frame(name = nodes), directed = FALSE)
g_ns <- make_g(e_ns, all_nodes); g_ss <- make_g(e_ss, all_nodes)

net_metrics <- function(g) {
  N <- vcount(g); E <- ecount(g)
  c(Edges = E, `Mean degree` = 2 * E / N,
    Density = edge_density(g),
    Modularity = modularity(cluster_fast_greedy(g)))
}
m_ns <- net_metrics(g_ns); m_ss <- net_metrics(g_ss)

# ============================================================
# F3b : dumbbell / arrow plot  (NS = 100% baseline)
# ============================================================
mlev <- c('Edges','Mean degree','Density','Modularity')
df <- data.frame(Metric = factor(mlev, levels = rev(mlev)),
                 NS = as.numeric(m_ns[mlev]),
                 SS = as.numeric(m_ss[mlev]))
df$rel_NS <- 100
df$rel_SS <- df$SS / df$NS * 100
df$dir    <- ifelse(df$SS < df$NS, 'Decrease', 'Increase')

long <- data.frame(
  Metric = rep(df$Metric, 2),
  Group  = factor(rep(c('NS','SS'), each = nrow(df)), levels = c('NS','SS')),
  rel    = c(df$rel_NS, df$rel_SS),
  val    = c(df$NS, df$SS))

fmt <- function(v) ifelse(v >= 10, sprintf('%.0f', v), sprintf('%.2f', v))

xmax <- max(long$rel) * 1.12
p_db <- ggplot() +
  geom_vline(xintercept = 100, linetype = 'dashed',
             color = 'gray60', linewidth = 0.4) +
  geom_segment(data = df,
               aes(y = Metric, yend = Metric, x = rel_NS, xend = rel_SS),
               color = 'gray55', linewidth = 1.2,
               arrow = arrow(length = unit(0.16, 'cm'), type = 'closed')) +
  geom_point(data = long, aes(y = Metric, x = rel, color = Group), size = 5.5) +
  geom_text(data = long,
            aes(y = Metric, x = rel, label = fmt(val), color = Group),
            vjust = -1.1, size = 3.1, fontface = 'bold', show.legend = FALSE) +
  annotate('text', x = 100, y = nrow(df) + 0.72, label = 'NS baseline',
           color = 'black', size = 3, fontface = 'italic', vjust = 0) +
  scale_color_manual(values = my_colors, name = NULL) +
  scale_x_continuous('Network metric relative to NS  (NS = 100%)',
                     limits = c(0, xmax),
                     breaks = c(0, 50, 100, 150),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(y = NULL,
       title = 'Network topology collapses from NS to SS') +
  coord_cartesian(ylim = c(0.5, nrow(df) + 0.5), expand = FALSE, clip = 'off') +
  theme_bw_pub(11) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = 'top',
        plot.margin = margin(24, 14, 6, 6),
        axis.text.y = element_text(face = 'bold', size = 11))

save_fig(p_db, 'F3b_network_metrics_dumbbell', OUT, 7.5, 3.8)
cat('Saved F3b dumbbell.\n'); print(df)
