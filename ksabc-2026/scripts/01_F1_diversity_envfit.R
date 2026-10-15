source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(vegan)
  library(reshape2)
  library(ggrepel)
})

DATA <- 'data'
OUT  <- 'output'

prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
mat <- t(as.matrix(gen[, samples]))
rownames(mat) <- samples

meta <- meta[samples, ]
meta$Group <- factor(meta$Group, levels = GROUP_ORDER)
pal <- my_colors

# ===== Alpha diversity (ACE, Shannon, Simpson) =====
ace <- estimateR(mat)['S.ACE', ]
shan <- diversity(mat, index = 'shannon')
simp <- diversity(mat, index = 'simpson')
alpha_df <- data.frame(SampleID = samples, Group = meta$Group,
                       ACE = ace, Shannon = shan, Simpson = simp)
alpha_long <- melt(alpha_df, id.vars = c('SampleID','Group'),
                   variable.name = 'Index', value.name = 'Value')

# Wilcoxon per index
pvals <- sapply(c('ACE','Shannon','Simpson'), function(idx) {
  d <- alpha_df[, idx]
  wilcox.test(d ~ alpha_df$Group)$p.value
})
lab_df <- data.frame(Index = names(pvals),
                     label = sprintf('italic(p)*" = %.3f"', pvals),
                     y = tapply(alpha_long$Value, alpha_long$Index,
                                function(v) max(v, na.rm = TRUE)) * 1.07)
lab_df$Index <- factor(lab_df$Index, levels = c('ACE','Shannon','Simpson'))
alpha_long$Index <- factor(alpha_long$Index, levels = c('ACE','Shannon','Simpson'))

p_alpha <- ggplot(alpha_long, aes(x = Group, y = Value)) +
  add_box_layers('Group') +
  geom_text(data = lab_df, aes(x = 1.5, y = y, label = label),
           inherit.aes = FALSE, size = 3.2, parse = TRUE) +
  facet_wrap(~ Index, scales = 'free_y', nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) +
  labs(x = NULL, y = 'Index value') +
  theme_bw_pub(11) +
  theme(strip.background = element_rect(fill = 'white', color = NA),
        strip.text = element_text(face = 'bold'),
        legend.position = 'none')

# ===== Beta diversity (Bray-Curtis PCoA) + envfit =====
bc <- vegdist(mat, method = 'bray')
pcoa <- cmdscale(bc, k = 4, eig = TRUE)
eig <- pcoa$eig
var_explained <- eig / sum(eig[eig > 0]) * 100
scores_pc <- as.data.frame(pcoa$points[, 1:2])
colnames(scores_pc) <- c('PC1','PC2')
scores_pc$Group <- meta$Group
scores_pc$SampleID <- samples

# Envfit with clinical vars
clin_vars <- c('Age','Burnish','Moisture','Sebum','Desquamation','Temperature','Redness')
env <- as.data.frame(meta[, clin_vars])
set.seed(42)
ef <- envfit(pcoa$points[, 1:2] ~ ., data = env, permutations = 999, na.rm = TRUE)
vec_df <- as.data.frame(scores(ef, 'vectors'))
colnames(vec_df) <- c('PC1','PC2')
vec_df$var <- rownames(vec_df)
vec_df$r2 <- ef$vectors$r
vec_df$p  <- ef$vectors$pvals

# Scale arrows
scl <- 0.45 * max(abs(scores_pc$PC1), abs(scores_pc$PC2)) /
       max(abs(c(vec_df$PC1, vec_df$PC2)))
vec_df$PC1s <- vec_df$PC1 * scl
vec_df$PC2s <- vec_df$PC2 * scl

# PERMANOVA
set.seed(42)
pn <- adonis2(bc ~ Group, data = meta, permutations = 999)
permR2 <- round(pn$R2[1], 3)
permP  <- pn$`Pr(>F)`[1]

# ANOSIM
set.seed(42)
ans <- anosim(bc, meta$Group, permutations = 999)
anR <- round(ans$statistic, 3)
anP <- ans$signif

stat_label <- bquote("PERMANOVA "*italic(R)^2*" = "*.(permR2)*", "*italic(p)*" = "*.(format(permP, scientific = FALSE)))

# Spider lines: centroid to each point per group
centroids <- aggregate(cbind(PC1, PC2) ~ Group, data = scores_pc, FUN = mean)
spider_df <- merge(scores_pc, centroids, by = 'Group', suffixes = c('','_c'))

# Label placement: significant labels stay at arrow tip; non-significant ones are
# nudged radially out to the plot periphery (with connector segments) to declutter.
ang   <- atan2(vec_df$PC2s, vec_df$PC1s)
rng_x <- max(abs(scores_pc$PC1)); rng_y <- max(abs(scores_pc$PC2))
# Push every label out along its arrow direction (significant a bit in, non-sig to
# the periphery), with connector segments back to the tip. The sample points are
# added as zero-label rows so ggrepel also avoids placing labels on top of points.
rad  <- ifelse(vec_df$p < 0.05, 0.85, 1.02)
vlab <- ifelse(vec_df$p < 0.05, sprintf('%s (r²=%.2f)', vec_df$var, vec_df$r2), vec_df$var)
repel_df <- rbind(
  data.frame(x = vec_df$PC1s, y = vec_df$PC2s, lab = vlab,
             nx = cos(ang) * rng_x * rad - vec_df$PC1s,
             ny = sin(ang) * rng_y * rad - vec_df$PC2s,
             fface = ifelse(vec_df$p < 0.05, 'bold', 'plain'),
             col   = ifelse(vec_df$p < 0.05, 'black', 'gray30')),
  data.frame(x = scores_pc$PC1, y = scores_pc$PC2, lab = '',
             nx = 0, ny = 0, fface = 'plain', col = 'black'))

p_pcoa <- ggplot(scores_pc, aes(PC1, PC2, color = Group, fill = Group)) +
  geom_segment(data = spider_df,
               aes(x = PC1_c, y = PC2_c, xend = PC1, yend = PC2,
                   color = Group),
               inherit.aes = FALSE, linewidth = 0.4, alpha = 0.4) +
  stat_ellipse(geom = 'polygon', alpha = 0.12, color = NA, level = 0.95) +
  geom_point(size = 3.2, alpha = 0.85,
             shape = 21, color = 'black', stroke = 0.4) +
  geom_point(data = centroids, aes(PC1, PC2, fill = Group),
             shape = 23, size = 4.5, color = 'black', stroke = 0.6,
             inherit.aes = FALSE) +
  geom_segment(data = vec_df,
               aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, 'cm')),
               color = ifelse(vec_df$p < 0.05, 'black', 'gray60'),
               linewidth = 0.6) +
  geom_text_repel(data = repel_df,
            aes(x = x, y = y, label = lab),
            inherit.aes = FALSE,
            nudge_x = repel_df$nx, nudge_y = repel_df$ny,
            fontface = repel_df$fface, color = repel_df$col,
            size = 3.1, seed = 7,
            box.padding = 0.5, point.padding = 0.3,
            force = 5, force_pull = 0.15, max.overlaps = Inf,
            min.segment.length = 0, segment.color = 'gray60', segment.size = 0.3) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(expand = expansion(mult = 0.13)) +
  scale_y_continuous(expand = expansion(mult = 0.13)) +
  labs(x = sprintf('PCo1 (%.1f%%)', var_explained[1]),
       y = sprintf('PCo2 (%.1f%%)', var_explained[2]),
       title = 'PCoA (Bray-Curtis) + envfit',
       subtitle = stat_label) +
  theme_bw_pub(11) +
  theme(panel.border = element_rect(color = 'grey20', fill = NA, linewidth = 0.9))

save_fig(p_alpha, 'F1a_alpha', OUT, 7.5, 3.2)
save_fig(p_pcoa,  'F1b_pcoa_envfit', OUT, 7.2, 5.8)

# Combined panel using gridExtra-free approach (base grid)
library(grid)
pdf(file.path(OUT, 'F1_diversity_combined.pdf'), width = 12, height = 4.2)
grid.newpage()
# alpha on left 55%, pcoa on right 45%
print(p_alpha, vp = viewport(x = 0.275, y = 0.5, width = 0.55, height = 1))
print(p_pcoa,  vp = viewport(x = 0.775, y = 0.5, width = 0.45, height = 1))
dev.off()
cat('Saved F1.\n')
cat('PERMANOVA R^2 =', permR2, ' p =', permP, '\n')
cat('ANOSIM R =', anR, ' p =', anP, '\n')
cat('Envfit significant (p<0.05):\n')
print(vec_df[vec_df$p < 0.05, c('var','r2','p')])
