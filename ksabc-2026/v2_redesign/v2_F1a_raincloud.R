source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(vegan)
  library(reshape2)
})

DATA <- 'data'
OUT  <- 'output/v2'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Self-contained half (flat) violin geom  — no extra packages
# ============================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
  setup_data = function(data, params) {
    # per-row only; do NOT split by group (would merge across facet panels)
    data$width <- data$width %||% params$width %||% 0.9
    data$xmin  <- data$x
    data$xmax  <- data$x + data$width / 2
    data
  },
  draw_group = function(data, panel_scales, coord) {
    data  <- transform(data, xminv = x, xmaxv = x + violinwidth * (xmax - x))
    left  <- transform(data, x = xminv); left  <- left[order(left$y), ]
    right <- transform(data, x = xmaxv); right <- right[order(-right$y), ]
    newdata <- rbind(left, right)
    newdata <- rbind(newdata, newdata[1, ])
    ggplot2:::ggname("geom_flat_violin",
      ggplot2::GeomPolygon$draw_panel(newdata, panel_scales, coord))
  },
  draw_key = draw_key_polygon,
  default_aes = aes(weight = 1, colour = "grey20", fill = "white",
                    linewidth = 0.5, alpha = NA, linetype = "solid"),
  required_aes = c("x", "y")
)
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
    position = "dodge", trim = TRUE, scale = "area", show.legend = NA,
    inherit.aes = TRUE, ...) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomFlatViolin,
        position = position, show.legend = show.legend,
        inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, ...))
}

# ============================================================
# Data + alpha diversity (identical computation to 01_F1)
# ============================================================
prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
mat <- t(as.matrix(gen[, samples])); rownames(mat) <- samples
meta <- meta[samples, ]
meta$Group <- factor(meta$Group, levels = GROUP_ORDER)

ace  <- estimateR(mat)['S.ACE', ]
shan <- diversity(mat, index = 'shannon')
simp <- diversity(mat, index = 'simpson')
alpha_df <- data.frame(SampleID = samples, Group = meta$Group,
                       ACE = ace, Shannon = shan, Simpson = simp)
alpha_long <- melt(alpha_df, id.vars = c('SampleID','Group'),
                   variable.name = 'Index', value.name = 'Value')

pvals <- sapply(c('ACE','Shannon','Simpson'), function(idx)
  wilcox.test(alpha_df[, idx] ~ alpha_df$Group)$p.value)
lab_df <- data.frame(Index = names(pvals),
                     label = sprintf('italic(p)*" = %.3f"', pvals),
                     y = tapply(alpha_long$Value, alpha_long$Index,
                                function(v) max(v, na.rm = TRUE)) * 1.10)
lev <- c('ACE','Shannon','Simpson')
lab_df$Index     <- factor(lab_df$Index, levels = lev)
alpha_long$Index <- factor(alpha_long$Index, levels = lev)

# ============================================================
# Raincloud:  half-violin (right) + box (right) + raw points (left)
# ============================================================
p_alpha <- ggplot(alpha_long, aes(x = Group, y = Value,
                                  fill = Group, color = Group)) +
  geom_flat_violin(position = position_nudge(x = 0.16), trim = TRUE,
                   width = 0.6, alpha = 0.55, color = NA) +
  geom_point(position = position_jitter(width = 0.055, seed = 1),
             size = 2, alpha = 0.6, shape = 16,
             aes(x = as.numeric(Group) - 0.12)) +
  geom_boxplot(position = position_nudge(x = 0.16), width = 0.10,
               outlier.shape = NA, alpha = 0.9,
               color = 'gray20', linewidth = 0.4) +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 3,
               color = 'black',
               position = position_nudge(x = 0.16)) +
  geom_text(data = lab_df, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, parse = TRUE) +
  facet_wrap(~ Index, scales = 'free_y', nrow = 1) +
  scale_fill_manual(values = my_colors, name = NULL) +
  scale_color_manual(values = my_colors, name = NULL, guide = 'none') +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.65))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.15))) +
  labs(x = NULL, y = 'Index value') +
  theme_bw_pub(11) +
  theme(strip.background = element_rect(fill = 'white', color = NA),
        strip.text = element_text(face = 'bold'),
        legend.position = 'none',
        panel.grid.major.x = element_blank())

save_fig(p_alpha, 'F1a_alpha_raincloud', OUT, 7.8, 3.3)
cat('Saved F1a raincloud.\n'); print(lab_df)
