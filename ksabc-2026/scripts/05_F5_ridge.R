source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(ggridges)
  library(dplyr)
  library(tidyr)
})

DATA <- 'data'
OUT  <- 'output'

prep <- readRDS(file.path(DATA, 'prep.rds'))
meta <- prep$meta
gen  <- prep$genus_rare

samples <- intersect(meta$SampleID, setdiff(colnames(gen), 'Genus'))
ra <- gen
ra[, samples] <- sweep(ra[, samples], 2, colSums(ra[, samples]), '/') * 100
rownames(ra) <- ra$Genus

# Pick a small set of key genera (driven by F4 results + biological story)
key_genera <- c('Lawsonella','Staphylococcus','Cutibacterium','Streptococcus',
                'Micrococcus','Latilactobacillus','Kocuria','Deinococcus')
key_genera <- key_genera[key_genera %in% rownames(ra)]

# Long format: SampleID, Genus, Sebum, RA
sebum <- meta[samples, 'Sebum']
df <- do.call(rbind, lapply(key_genera, function(g) {
  data.frame(SampleID = samples,
             Genus    = g,
             Group    = meta[samples, 'Group'],
             Sebum    = sebum,
             RA       = as.numeric(ra[g, samples]))
}))

# Sebum quartile (across all samples — NS very low, SS very high; use 4 bins)
brks <- quantile(df$Sebum, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
df$SebumQ <- cut(df$Sebum, breaks = brks, include.lowest = TRUE,
                 labels = c('Q1 (lowest)','Q2','Q3','Q4 (highest)'))

# Order genera by Spearman with Sebum (descending)
ordr <- df %>%
  group_by(Genus) %>%
  summarise(r = suppressWarnings(cor(Sebum, RA, method = 'spearman')),
            mean_RA = mean(RA)) %>%
  arrange(desc(r))
df$Genus <- factor(df$Genus, levels = rev(ordr$Genus))

pal_q <- c('Q1 (lowest)'  = '#002060',
           'Q2'           = '#0070C0',
           'Q3'           = '#E79453',
           'Q4 (highest)' = '#C0000D')

# Use log10(RA + 0.01) for visual spread
df$logRA <- log10(df$RA + 0.01)

p_ridge <- ggplot(df, aes(x = logRA, y = Genus, fill = SebumQ)) +
  geom_density_ridges(alpha = 0.75, scale = 1.2, color = 'gray20',
                      linewidth = 0.25, bandwidth = 0.35) +
  scale_fill_manual(values = pal_q, name = 'Sebum quartile') +
  scale_x_continuous(breaks = log10(c(0.01, 0.1, 1, 10, 100) + 0.01),
                     labels = c('0','0.1','1','10','100')) +
  labs(x = 'Relative abundance (%)',
       y = NULL,
       title = 'Genus abundance distributions across sebum quartiles') +
  theme_ridges(grid = FALSE) +
  theme(axis.text = element_text(color = 'black'),
        axis.text.y = element_text(face = 'italic', color = 'black'),
        legend.position = 'right',
        plot.title = element_text(face = 'bold', size = 13),
        plot.subtitle = element_text(size = 9, color = 'black'))

# Companion: Lawsonella vs Sebum scatter with loess fit
law  <- df[df$Genus == 'Lawsonella', ]
rval <- round(cor(law$Sebum, law$RA, method = 'spearman'), 2)
p_scatter <- ggplot(law, aes(Sebum, RA)) +
  geom_smooth(method = 'loess', se = TRUE, color = '#C0000D',
              fill = '#C0000D', alpha = 0.18, linewidth = 0.8) +
  geom_point(aes(fill = Group), shape = 21, color = 'black',
             stroke = 0.4, size = 3, alpha = 0.9) +
  scale_fill_manual(values = my_colors) +
  labs(x = 'Sebum (clinical)',
       y = expression(paste(italic("Lawsonella"), " RA (%)")),
       title = expression(paste(bolditalic("Lawsonella"), bold(" vs. sebum"))),
       subtitle = bquote("Spearman "*italic(r)*" = "*.(rval)*", "*italic(p)*" < 0.001")) +
  theme_bw_pub(11) +
  theme(plot.title = element_text(size = 13))

save_fig(p_ridge,   'F5a_ridge_sebumQ',      OUT, 8.5, 5.5)
save_fig(p_scatter, 'F5b_law_sebum_scatter', OUT, 5,   4)

# Combined panel
library(grid)
pdf(file.path(OUT, 'F5_combined.pdf'), width = 13.5, height = 5.5)
grid.newpage()
print(p_ridge,   vp = viewport(x = 0.36, y = 0.5, width = 0.66, height = 1))
print(p_scatter, vp = viewport(x = 0.83, y = 0.5, width = 0.34, height = 1))
dev.off()

cat('Saved F5.\n')
print(ordr)
