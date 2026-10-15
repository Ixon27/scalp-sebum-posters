#!/usr/bin/env python
"""
15  Redness(홍반) 가 다양성 축의 주인인가 — 수분과 갈라낸다.

14 에서 나온 것 (나이 통제, BH-FDR):
    Shannon  ~ Redness +0.48*  vs  Moisture +0.47*
    Observed ~ Redness +0.57*  vs  Moisture +0.23 (비유의)
  Redness 와 Moisture 는 서로 상관이 없다(+0.123, 비유의). 즉 별개의 변수 둘이 같은 것을 설명한다.

가릴 것
  (1) 서로를 통제하면 누가 남는가 (나이 + 상대변수 동시 통제)
  (2) 조성 전체(Bray-Curtis)의 고유 설명력은 어떤가 — 두 축 분산분할에 Redness 를 넣으면
  (3) 이상치에 끌려간 것은 아닌가 — SS01 의 Redness 34.0 은 나머지(92~220)에서 크게 벗어난다
  (4) Redness 는 군(NS/SS)을 가르지 못한다(14-D: q=0.43). 그렇다면 무엇을 재는 값인가

결과가 어떻게 나오든 그대로 보고한다. Redness 가 이기면 초록의 'hydration governs diversity' 를 고쳐야 한다.
"""
import csv, collections
import numpy as np
from scipy import stats
from scipy.spatial.distance import braycurtis

SEED, NPERM = 20260825, 9999
rng = np.random.default_rng(SEED)
ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
rep = []
def W(*a):
    s = "".join(str(x) for x in a); rep.append(s); print(s, flush=True)

md = {}
with open(f"{ROOT}/output/01_metadata.tsv") as f:
    r = csv.DictReader(f, delimiter="\t"); idc = r.fieldnames[0]
    for row in r: md[row[idc]] = row
gc = collections.OrderedDict()
with open(f"{ROOT}/output/01_genus_counts.tsv") as f:
    r = csv.reader(f, delimiter="\t"); samples = next(r)[1:]
    for row in r: gc[row[0]] = np.array([float(x) for x in row[1:]])
genera = list(gc); gi = {g: i for i, g in enumerate(genera)}
G = np.vstack([gc[k] for k in genera])
IDX = ['Sebum', 'Moisture', 'Age', 'Burnish', 'Desquamation', 'Temperature', 'Redness']
X = {k: np.array([float(md[s][k]) for s in samples]) for k in IDX}
grp = np.array([md[s]['Group'] for s in samples])


def shannon(v):
    p = v[v > 0] / v.sum(); return float(-(p * np.log(p)).sum())
def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d: return v.astype(float)
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False), minlength=len(v)).astype(float)

depth = int(G.sum(axis=0).min())
R = np.column_stack([rarefy(G[:, j], depth) for j in range(G.shape[1])])
rel = R / R.sum(axis=0)
sh = np.array([shannon(R[:, j]) for j in range(R.shape[1])])
obs_n = (R > 0).sum(axis=0).astype(float)
dom = rel.max(axis=0)
lr_SL = np.log((R[gi['Staphylococcus']] + .5) / (R[gi['Lawsonella']] + .5))

W("=" * 78)
W("15  Redness 가 다양성 축의 주인인가 — 수분과 갈라낸다")
W(f"    시료 {len(samples)}개 · 깊이 {depth:,} · 순열 {NPERM}회 · 난수 {SEED}")
W("=" * 78)

W("\n### 0. Redness 값의 분포부터")
rd = X['Redness']
W(f"  범위 {rd.min():.2f} ~ {rd.max():.2f} · 중앙 {np.median(rd):.2f} · 평균 {rd.mean():.2f}")
lo = [(samples[i], rd[i], grp[i]) for i in np.argsort(rd)[:3]]
hi = [(samples[i], rd[i], grp[i]) for i in np.argsort(rd)[-3:]]
W(f"  최저 3: {lo}")
W(f"  최고 3: {hi}")
z = (rd - rd.mean()) / rd.std()
out = [(samples[i], rd[i], z[i]) for i in range(len(rd)) if abs(z[i]) > 2]
W(f"  |z|>2 이상치: {out if out else '없음'}")


def pcor_multi(x, y, Zlist, n=NPERM):
    """Zlist 전부 통제한 순위 편상관 + 순열"""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    Z = np.c_[np.ones(len(x))] if not Zlist else np.c_[np.ones(len(x)), *[stats.rankdata(z) for z in Zlist]]
    f = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    a, b = f(rx), f(ry)
    o = np.corrcoef(a, b)[0, 1]
    null = np.array([np.corrcoef(rng.permutation(a), b)[0, 1] for _ in range(n)])
    return o, (np.sum(np.abs(null) >= abs(o)) + 1) / (n + 1)


W("\n### 1. 서로를 통제하면 누가 남는가")
TARG = [("Shannon", sh), ("Observed", obs_n), ("우점도", dom)]
for tn, y in TARG:
    W(f"\n  ── {tn}")
    r1, p1 = pcor_multi(X['Moisture'], y, [X['Age']])
    r2, p2 = pcor_multi(X['Redness'], y, [X['Age']])
    W(f"     {'수분 (나이만 통제)':<32}rho={r1:+.3f}  p={p1:.4f}")
    W(f"     {'홍반 (나이만 통제)':<32}rho={r2:+.3f}  p={p2:.4f}")
    r3, p3 = pcor_multi(X['Moisture'], y, [X['Age'], X['Redness']])
    r4, p4 = pcor_multi(X['Redness'], y, [X['Age'], X['Moisture']])
    W(f"     {'수분 (나이+홍반 통제)':<32}rho={r3:+.3f}  p={p3:.4f}   {'살아남음' if p3<0.05 else '탈락'}")
    W(f"     {'홍반 (나이+수분 통제)':<32}rho={r4:+.3f}  p={p4:.4f}   {'살아남음' if p4<0.05 else '탈락'}")

W("\n### 2. 이상치 SS01(Redness 34.0) 을 빼면")
keep = np.array([s != 'SS01' for s in samples])
for tn, y in TARG:
    r2, p2 = pcor_multi(X['Redness'][keep], y[keep], [X['Age'][keep]])
    r1, p1 = pcor_multi(X['Moisture'][keep], y[keep], [X['Age'][keep]])
    W(f"  {tn:<10} 홍반 rho={r2:+.3f} p={p2:.4f}    수분 rho={r1:+.3f} p={p1:.4f}   (n={keep.sum()})")

W("\n### 3. 피지 축(로그비)에는 홍반이 붙나 — 두 축 구조가 깨지는지")
for nm, z in (("나이만 통제", [X['Age']]), ("나이+피지 통제", [X['Age'], X['Sebum']])):
    r, p = pcor_multi(X['Redness'], lr_SL, z)
    W(f"  홍반 ~ log(Staph/Laws) ({nm:<14}) rho={r:+.3f}  p={p:.4f}")
r, p = pcor_multi(X['Sebum'], lr_SL, [X['Age'], X['Redness']])
W(f"  피지 ~ log(Staph/Laws) (나이+홍반 통제) rho={r:+.3f}  p={p:.4f}   {'유지' if p<0.05 else '탈락'}")

W("\n### 4. 조성 전체(Bray-Curtis) 고유 설명력 — 홍반을 넣으면")
n = len(samples)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        D[i, j] = D[j, i] = braycurtis(rel[:, i], rel[:, j])
A = -0.5 * D ** 2; Cn = np.eye(n) - np.ones((n, n)) / n
Gm = Cn @ A @ Cn; tot = np.trace(Gm)


def r2_partial(target, covars):
    """covars 를 통제한 뒤 target 의 고유 R²"""
    def rk(v):
        v = stats.rankdata(v).astype(float); return (v - v.mean()) / v.std()
    t = rk(target)
    if covars:
        Z = np.c_[np.ones(n), *[rk(c) for c in covars]]
        t = t - Z @ np.linalg.lstsq(Z, t, rcond=None)[0]
        if t.std() == 0: return np.nan
        t = t / t.std()
    H = np.outer(t, t) / (t @ t)
    return np.trace(H @ Gm @ H) / tot


combos = [("피지", X['Sebum'], [X['Moisture'], X['Age'], X['Redness']]),
          ("수분", X['Moisture'], [X['Sebum'], X['Age'], X['Redness']]),
          ("나이", X['Age'], [X['Sebum'], X['Moisture'], X['Redness']]),
          ("홍반", X['Redness'], [X['Sebum'], X['Moisture'], X['Age']])]
W(f"  {'변수':<8}{'고유 R²':>10}{'p':>9}   (나머지 셋 전부 통제)")
for nm, t, cv in combos:
    o = r2_partial(t, cv)
    null = np.array([r2_partial(rng.permutation(t), cv) for _ in range(999)])
    p = (np.sum(null >= o) + 1) / 1000
    W(f"  {nm:<8}{o:>10.4f}{p:>9.4f}" + ("   **유의**" if p < 0.05 else ""))

W("\n### 5. 홍반은 무엇을 재는 값인가 — 참고")
W(f"  NS 중앙 {np.median(rd[grp=='NS']):.1f} vs SS 중앙 {np.median(rd[grp=='SS']):.1f}"
  f"  (Mann-Whitney p={stats.mannwhitneyu(rd[grp=='NS'], rd[grp=='SS']).pvalue:.4f}, 14-D 에서 q=0.43 비유의)")
for k in ['Sebum', 'Moisture', 'Age', 'Temperature', 'Desquamation', 'Burnish']:
    r, p = pcor_multi(rd, X[k], [])
    W(f"  홍반 ~ {k:<14} rho={r:+.3f}  p={p:.4f}" + ("  *" if p < 0.05 else ""))

W("\n" + "=" * 78)
open(f"{ROOT}/output/15_redness.txt", "w").write("\n".join(rep))
print(f"\n저장: {ROOT}/output/15_redness.txt")
