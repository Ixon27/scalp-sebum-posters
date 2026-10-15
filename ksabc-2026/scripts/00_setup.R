# Common setup for KSABC figures
# Adheres to user's microbiome 분석 visual style guide
.libPaths(c('~/R_libs', .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
})

# ===== Palettes (from user's microbiome_분석 standard) =====
my_colors <- c('NS' = '#002060', 'SS' = '#C0000D')

bar_palette <- c('#002060','#0070C0','#C00000','#00B0F0','#E4C868',
                 '#BF9000','#00B050','#92D050','#E79453','#C84E40',
                 'grey85')  # last = Others

# Divergent palette for heatmaps (blue - yellow - red)
heatmap_div <- c('#2b83ba','#ffffbf','#d7191c')

# ===== Themes (per user's theme_pub / theme_bw_pub) =====
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = 'sans'),
      axis.title = element_text(face = 'bold'),
      axis.text = element_text(color = 'black'),
      plot.title = element_text(face = 'bold', size = base_size + 2),
      plot.subtitle = element_text(size = base_size - 1.5, color = 'black'),
      legend.title = element_text(face = 'bold'),
      legend.background = element_blank()
    )
}

theme_bw_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = 'sans'),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = 'bold'),
      axis.text = element_text(color = 'black'),
      plot.title = element_text(face = 'bold', size = base_size + 2),
      plot.subtitle = element_text(size = base_size - 1.5, color = 'black'),
      legend.title = element_text(face = 'bold')
    )
}

# ===== Helpers =====
# 2-group boxplot block (boxplot + jitter + mean diamond, fills via my_colors)
add_box_layers <- function(group_col = 'Group') {
  list(
    geom_boxplot(aes(fill = !!sym(group_col)),
                 alpha = 0.8, outlier.shape = NA, width = 0.5,
                 color = 'gray15', linewidth = 0.4),
    geom_jitter(aes(color = !!sym(group_col)),
                width = 0.2, size = 2.5, alpha = 0.5),
    stat_summary(fun = mean, geom = 'point', shape = 18,
                 size = 4, color = 'black'),
    scale_fill_manual(values = my_colors, name = NULL),
    scale_color_manual(values = my_colors, name = NULL, guide = 'none')
  )
}
# rlang :: sym for !!sym usage
if (!requireNamespace('rlang', quietly = TRUE)) {
  sym <- function(x) as.name(x)
} else {
  sym <- rlang::sym
}

# Wilcoxon bracket between two x positions on a ggplot
add_bracket <- function(p, y_top, p_val, x1 = 1, x2 = 2, tick_frac = 0.015,
                        label_offset = 0.04, lwd = 0.6) {
  ylim <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
  if (is.null(ylim)) ylim <- c(0, y_top)
  y_range <- diff(ylim)
  tick <- y_range * tick_frac
  lab <- if (p_val < 0.001) '***' else if (p_val < 0.01) '**' else
         if (p_val < 0.05)  '*'   else 'ns'
  lab_full <- sprintf('%s  p = %.3g', lab, p_val)
  p +
    annotate('segment', x = x1, xend = x2, y = y_top, yend = y_top, linewidth = lwd) +
    annotate('segment', x = x1, xend = x1, y = y_top, yend = y_top - tick, linewidth = lwd) +
    annotate('segment', x = x2, xend = x2, y = y_top, yend = y_top - tick, linewidth = lwd) +
    annotate('text', x = (x1 + x2) / 2, y = y_top + y_range * label_offset,
             label = lab_full, fontface = 'bold', size = 3.4)
}

# Save fig in both PDF (cairo, fonttype-safe) and PNG (cairo)
save_fig <- function(plot, name, out_dir, w, h, dpi = 300) {
  ggsave(file.path(out_dir, paste0(name, '.pdf')),
         plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, '.png')),
         plot, width = w, height = h, dpi = dpi, type = 'cairo')
}

GROUP_ORDER <- c('NS','SS')

# 결과 폴더
dir.create('output', showWarnings = FALSE)
