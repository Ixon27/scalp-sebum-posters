#!/usr/bin/env python
"""
18  신뢰구간 포함률 시뮬레이션 — 17 에서 쓴 구간 공식이 우리 쓰임새에서도 95% 를 지키는가.

왜: Bonett & Wright(2000) 의 Spearman 분산 (1 + r²/2)/(n − 3) 은 **단순 상관**용이다
    (2차 출처 2건으로 공식 확인: Ruscio 2008 식 7, Bishara & Hittner 2017 식 16).
    17 은 이를 **편상관**으로 확장해 자유도를 n − 3 − k 로 줄였고(k = 통제 변수 수),
    군(0/1)이 x 를 완전히 가르는 경우(피지)까지 썼다. 이 확장은 원문·2차 출처 범위 밖이라 직접 확인한다.

설계
  참값 = 같은 추정식(순위 → 통제변수 순위·군에 회귀 → 잔차 상관)을 N = 400,000 에 적용한 값
  표본 n ∈ {15, 27}, 오차 상관 ρe ∈ {0, 0.3, 0.6}, 반복 4,000회, 난수 20260825
  시나리오
    k=0  x, y 이변량 정규
    k=1  나이 같은 연속 교란 z 가 x·y 둘 다에 들어간다 (x = 0.6z + e, y = 0.6z + e)
    k=2  z + 군 g(0/1, 비율 11/27). x = 3.0g + 0.4z + e  (군이 x 를 거의 완전히 가른다 — 피지처럼)
                                     y = 0.8g + 0.4z + e
  비교하는 구간
    BW-z   : z ± 1.96 · sqrt((1 + r²/2)/(n − 3 − k))     ← 17 이 쓴 것
    BW-t   : 같은 SE, 임계값 t(0.975, n − 3 − k)          ← Ruscio 2008 이 언급한 변형
    Fisher : z ± 1.96 / sqrt(n − 3 − k)                    ← Pearson 공식
  순위는 단조 변환에 불변이라 정규 주변분포 하나로 모든 단조 모양을 대표한다(동점은 제외).

실행: python3 analysis/17_small_n/18_ci_coverage.py   (로컬 수십 초, 병렬 없음)
"""
import numpy as np
from scipy import stats

ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
RNG = np.random.default_rng(20260825)
REPS, NBIG = 4000, 400_000
LOG = []


def say(s=""):
    print(s)
    LOG.append(s)


def pcor(x, y, covs):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    Z = np.column_stack([np.ones(len(rx))] + [c if set(np.unique(c)) <= {0.0, 1.0}
                                              else stats.rankdata(c) for c in covs])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def draw(n, k, rho_e):
    e = RNG.multivariate_normal([0, 0], [[1, rho_e], [rho_e, 1]], size=n)
    if k == 0:
        return e[:, 0], e[:, 1], []
    z = RNG.normal(size=n)
    if k == 1:
        return 0.6 * z + e[:, 0], 0.6 * z + e[:, 1], [z]
    g = (RNG.random(n) < 11 / 27).astype(float)
    if g.sum() < 2 or g.sum() > n - 2:      # 한 군이 거의 비면 다시 뽑는다
        return draw(n, k, rho_e)
    return 3.0 * g + 0.4 * z + e[:, 0], 0.8 * g + 0.4 * z + e[:, 1], [z, g]


say("=" * 78)
say("18  신뢰구간 포함률 — 편 Spearman + Bonett–Wright 분산 (n − 3 − k)")
say(f"    반복 {REPS}회 · 참값 N = {NBIG:,} · 난수 20260825 · 명목 포함률 0.95")
say("    포함률 몬테카를로 오차: ±0.007 (95%, 4,000회)")
say("=" * 78)
say(f"\n  {'k':>2} {'n':>3} {'ρe':>4} {'참값':>7} | {'BW-z':>6} {'BW-t':>6} {'Fisher':>6} | {'BW-z 평균폭':>9}")
rows = []
for k in (0, 1, 2):
    for rho_e in (0.0, 0.3, 0.6):
        xb, yb, cb = draw(NBIG, k, rho_e)
        truth = pcor(xb, yb, cb)
        for n in (15, 27):
            df = n - 3 - k
            tcrit = stats.t.ppf(0.975, df)
            hit = np.zeros(3); width = 0.0
            for _ in range(REPS):
                x, y, c = draw(n, k, rho_e)
                r = np.clip(pcor(x, y, c), -0.9999, 0.9999)
                z = np.arctanh(r)
                se_bw = np.sqrt((1 + r * r / 2) / df)
                se_f = 1 / np.sqrt(df)
                lo1, hi1 = np.tanh(z - 1.96 * se_bw), np.tanh(z + 1.96 * se_bw)
                lo2, hi2 = np.tanh(z - tcrit * se_bw), np.tanh(z + tcrit * se_bw)
                lo3, hi3 = np.tanh(z - 1.96 * se_f), np.tanh(z + 1.96 * se_f)
                hit += [lo1 <= truth <= hi1, lo2 <= truth <= hi2, lo3 <= truth <= hi3]
                width += hi1 - lo1
            cov = hit / REPS
            rows.append((k, n, rho_e, truth, *cov, width / REPS))
            say(f"  {k:>2} {n:>3} {rho_e:>4.1f} {truth:>+7.3f} | {cov[0]:>6.3f} {cov[1]:>6.3f} "
                f"{cov[2]:>6.3f} | {width / REPS:>9.3f}")

cov = np.array([r[4:7] for r in rows])
say(f"\n  전체 범위 — BW-z {cov[:, 0].min():.3f}–{cov[:, 0].max():.3f} · "
    f"BW-t {cov[:, 1].min():.3f}–{cov[:, 1].max():.3f} · Fisher {cov[:, 2].min():.3f}–{cov[:, 2].max():.3f}")
for k in (0, 1, 2):
    m = np.array([r[0] == k for r in rows])
    say(f"  k={k}: BW-z 평균 {cov[m, 0].mean():.3f} · BW-t {cov[m, 1].mean():.3f} · Fisher {cov[m, 2].mean():.3f}")

with open(f"{ROOT}/output/18_ci_coverage.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
print("\n저장: output/18_ci_coverage.txt")
