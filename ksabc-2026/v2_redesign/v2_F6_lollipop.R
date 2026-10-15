source('scripts/00_setup.R')   # ksabc-2026/ 에서 실행
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(reshape2)
})

clean_desc <- function(x) {
  x <- gsub('&beta;','β',x);  x <- gsub('&alpha;','α',x)
  x <- gsub('&gamma;','γ',x); x <- gsub('&delta;','δ',x)
  x <- gsub('&Delta;','Δ',x); x <- gsub('&epsilon;','ε',x)
  x <- gsub('&omega;','ω',x); x <- gsub('&mu;','µ',x)
  x <- gsub('&amp;','&',x);   x <- gsub('<[^>]+>','',x)
  x
}

DATA <- 'data'
OUT  <- 'output/v2'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P2 <- file.path(DATA, 'picrust2_out')

prep <- readRDS(file.path(DATA, 'prep.rds')); meta <- prep$meta

pa_file <- file.path(P2, 'pathways_out', 'path_abun_unstrat.tsv.gz')
stopifnot(file.exists(pa_file))
pa <- read.table(gzfile(pa_file), header = TRUE, sep = '\t',
                 quote = '', check.names = FALSE, comment.char = '')
colnames(pa)[1] <- 'pathway'

desc_files <- c(Sys.getenv('PICRUST2_METACYC_DESC')   # PICRUSt2 설치 폴더의 default_files/description_mapfiles/metacyc_pathways_info.txt.gz (없으면 설명 없이 진행))
desc_files <- desc_files[file.exists(desc_files)]
if (length(desc_files) > 0) {
  desc <- read.table(gzfile(desc_files[1]), sep = '\t', header = FALSE,
                     quote = '', comment.char = '',
                     col.names = c('pathway','description'))
  pa <- merge(pa, desc, by = 'pathway', all.x = TRUE)
} else pa$description <- pa$pathway
pa$description <- clean_desc(pa$description)

samples <- intersect(meta$SampleID, colnames(pa)); meta <- meta[samples, ]
mat <- as.matrix(pa[, samples]); keep <- rowSums(mat) > 0
mat <- mat[keep, ]; info <- pa[keep, c('pathway','description')]
ra <- sweep(mat, 2, colSums(mat), '/') * 100

grp <- meta$Group; ns_s <- samples[grp == 'NS']; ss_s <- samples[grp == 'SS']
stats <- data.frame(
  pathway = info$pathway, description = info$description,
  mean_NS = rowMeans(ra[, ns_s]), mean_SS = rowMeans(ra[, ss_s]),
  log2FC = log2((rowMeans(ra[, ss_s]) + 1e-4) / (rowMeans(ra[, ns_s]) + 1e-4)),
  p = apply(ra, 1, function(v) tryCatch(
    wilcox.test(v[ns_s], v[ss_s], exact = FALSE)$p.value,
    error = function(e) NA)))
stats$padj <- p.adjust(stats$p, method = 'BH')

lipid_terms <- c('lipid','steroid','sterol','cholesterol','sebum','fatty acid',
                 'fatty-acid','sphingo','glycerolipid','glycerophos','triacyl',
                 'wax','bile acid','phospholipid','arachidon','oleate',
                 'palmitoleate','vaccenate','stearate','palmitate','linoleate',
                 'isoprenoid','squalene')
hit <- Reduce(`|`, lapply(lipid_terms,
              function(t) grepl(t, stats$description, ignore.case = TRUE)))
lipid <- stats[hit, ]
lipid$score <- abs(lipid$log2FC) * pmax(-log10(lipid$p), 0)
lipid <- lipid[order(-lipid$score), ]
top <- head(lipid, 20); top <- top[order(top$log2FC), ]
top$short <- substr(top$description, 1, 70)
top$short <- ifelse(nchar(top$description) > 70, paste0(top$short, '...'), top$short)
top$short <- factor(top$short, levels = top$short)
top$direction <- ifelse(top$log2FC > 0, 'SS-enriched', 'NS-enriched')
top$sig <- ifelse(top$p < 0.05, '*', '')

# ============================================================
# F6 : diverging lollipop  (replaces grouped bars)
# ============================================================
# Italicise organism names (E. coli, Chlamydia) inside pathway labels via plotmath
italicize_org <- function(x) {
  x <- gsub('E. coli',   '"*italic("E. coli")*"',   x, fixed = TRUE)
  x <- gsub('Chlamydia', '"*italic("Chlamydia")*"', x, fixed = TRUE)
  paste0('"', x, '"')
}

p_lolli <- ggplot(top, aes(x = log2FC, y = short, color = direction)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = log2FC, y = short, yend = short),
               linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 3.6) +
  geom_text(aes(label = sprintf('italic(p)*" = %.3g %s"', p, sig)),
            parse = TRUE,
            hjust = ifelse(top$log2FC > 0, -0.18, 1.18),
            size = 2.9, color = 'black') +
  scale_color_manual(values = c('NS-enriched' = '#002060',
                                'SS-enriched' = '#C0000D'), name = 'Direction') +
  scale_x_continuous(expand = expansion(mult = c(0.28, 0.28))) +
  scale_y_discrete(labels = function(b) parse(text = italicize_org(b))) +
  labs(x = expression(log[2]~`(SS / NS)`), y = NULL,
       title = 'Lipid / sebum-related MetaCyc pathway differences (top 20)') +
  theme_bw_pub(11) +
  theme(panel.grid.major.y = element_line(color = 'gray92', linewidth = 0.3),
        legend.position = 'top',
        axis.text = element_text(color = 'black'))

save_fig(p_lolli, 'F6_lipid_pathways_lollipop', OUT, 11, 7)
cat('Saved F6 lollipop.\n')
print(head(top[, c('description','log2FC','p')], 10))
