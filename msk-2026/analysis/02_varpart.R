# S1-2  조성 전체에서 피지 축과 수분 축의 '고유' 설명력
#
# 옛 초록은 envfit 으로 변수를 하나씩 따로 재서 sebum r2=0.53 < moisture r2=0.56 이라는
# 이유로 "a leading driver" 로 후퇴했다. 변수를 따로 재면 공유 설명력이 양쪽에 중복 계상된다.
# 여기서는 분산분할(varpart)과 부분 제약서열분석(partial dbRDA)으로 고유분을 분리한다.
# 연속 변수를 그대로 쓴다 — 중앙값 이분은 정보를 버린다.

suppressPackageStartupMessages(library(vegan))
set.seed(20260825)

root <- if (dir.exists("output")) "." else ".."
out  <- file.path(root, "output")

g  <- read.delim(file.path(out, "01_genus_counts.tsv"), row.names = 1, check.names = FALSE)
md <- read.delim(file.path(out, "01_metadata.tsv"),     row.names = 1, check.names = FALSE)
g  <- g[, rownames(md), drop = FALSE]
comm <- t(g)                       # 행 = 시료
ra   <- comm / rowSums(comm)

con <- file(file.path(out, "02_varpart.txt"), open = "wt", encoding = "UTF-8")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }

w(strrep("=", 78))
w("S1-2  조성에서 두 축의 고유 설명력 (n=", nrow(comm), " 시료, ", ncol(comm), " 속)")
w(strrep("=", 78))

# ---- 거리 두 가지: Bray-Curtis(상대존재량) / Aitchison(조성 인지)
clr <- function(x, pseudo = 0.5) { y <- log(x + pseudo); sweep(y, 1, rowMeans(y), "-") }
D_bray <- vegdist(ra, method = "bray")
D_ait  <- dist(clr(comm))

for (nm in c("Bray-Curtis", "Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_bray else D_ait
  w("\n### ", nm)

  # 1) 각 변수 단독 (옛 envfit 방식에 대응)
  w("\n  [단독] 변수를 따로 잰 설명력 — 공유분이 중복 계상된다")
  for (v in c("Sebum", "Moisture", "Age")) {
    m <- dbrda(D ~ md[[v]]); a <- anova(m, permutations = 9999)
    w(sprintf("    %-9s R2 = %.4f   F = %5.2f   p = %.4f", v,
              RsquareAdj(m)$r.squared, a$F[1], a$`Pr(>F)`[1]))
  }

  # 2) 부분 dbRDA — 나머지를 통제한 고유 설명력
  w("\n  [고유] 나머지 두 변수를 통제한 뒤 남는 설명력 (partial dbRDA)")
  for (v in c("Sebum", "Moisture", "Age")) {
    others <- setdiff(c("Sebum", "Moisture", "Age"), v)
    f <- as.formula(sprintf("D ~ %s + Condition(%s)", v, paste(others, collapse = " + ")))
    m <- dbrda(f, data = md); a <- anova(m, permutations = 9999)
    w(sprintf("    %-9s R2(고유) = %.4f   F = %5.2f   p = %.4f", v,
              RsquareAdj(m)$r.squared, a$F[1], a$`Pr(>F)`[1]))
  }

  # 3) 분산분할
  vp <- varpart(D, ~Sebum, ~Moisture, ~Age, data = md)
  fr <- vp$part$indfract
  lbl <- c("Sebum 고유", "Moisture 고유", "Age 고유",
           "Sebum∩Moisture", "Moisture∩Age", "Sebum∩Age", "3자 공유", "설명 안 됨")
  w("\n  [분산분할] adjusted R2")
  for (i in seq_len(nrow(fr))) w(sprintf("    %-16s %+.4f", lbl[i], fr$Adj.R.square[i]))

  # 4) 군 효과 대조
  m <- dbrda(D ~ Group, data = md); a <- anova(m, permutations = 9999)
  w(sprintf("\n  [참고] Group(NS/SS) 단독  R2 = %.4f  p = %.4f", RsquareAdj(m)$r.squared, a$`Pr(>F)`[1]))
  m2 <- dbrda(D ~ Group + Condition(Moisture + Age), data = md); a2 <- anova(m2, permutations = 9999)
  w(sprintf("         Group | Moisture+Age 통제   R2 = %.4f  p = %.4f",
            RsquareAdj(m2)$r.squared, a2$`Pr(>F)`[1]))
}

# ---- 분산의 등질성 (PERMANOVA 해석의 전제)
w("\n### 분산 등질성 검정 (betadisper) — 군간 차이가 '퍼짐의 차이'는 아닌지")
for (nm in c("Bray-Curtis", "Aitchison")) {
  D <- if (nm == "Bray-Curtis") D_bray else D_ait
  bd <- betadisper(D, md$Group); pt <- permutest(bd, permutations = 9999)
  w(sprintf("    %-12s F = %5.2f  p = %.4f  (NS 평균거리 %.3f / SS %.3f)",
            nm, pt$tab$F[1], pt$tab$`Pr(>F)`[1],
            mean(bd$distances[md$Group == "NS"]), mean(bd$distances[md$Group == "SS"])))
}

close(con)
cat("\n저장: output/02_varpart.txt\n")
