#!/usr/bin/env python
"""
14  임상 지표 7개를 전부 검정한다 — 정말 두 축뿐인가.

발단: 초록에 "seven quantitative scalp indices" 라고 썼는데, 실제 분석에 쓴 것은
      Sebum · Moisture · Age 셋뿐이다. 나머지 넷(Burnish 윤기, Desquamation 각질,
      Temperature, Redness)은 한 번도 검정되지 않았다.

두 가지 결과 중 하나가 나온다.
  (a) 넷 중 하나가 피지·수분만큼 또는 그 이상 설명한다 -> 두 축 프레임이 불완전하다. 논지 수정.
  (b) 넷 다 비유의 -> "일곱 개를 전부 재봤고 둘만 구조를 만든다" 가 되어 주장이 강해진다.

미리 정한 규칙
  - 다중검정: 지표군(7개) 안에서 BH-FDR.
  - 나이는 통제 변수이자 검정 대상 둘 다로 다룬다(나이 자신을 볼 때는 통제하지 않는다).
  - 결과가 (a) 로 나와도 그대로 보고한다.
"""
import csv, itertools, collections
import numpy as np
from scipy import stats

SEED, NPERM = 20260825, 9999
rng = np.random.default_rng(SEED)
ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
rep = []
def W(*a):
    s = "".join(str(x) for x in a); rep.append(s); print(s, flush=True)


md = {}
with open(f"{ROOT}/output/01_metadata.tsv") as f:
    r = csv.DictReader(f, delimiter="\t"); idc = r.fieldnames[0]
    for row in r:
        md[row[idc]] = row

gc = collections.OrderedDict()
with open(f"{ROOT}/output/01_genus_counts.tsv") as f:
    r = csv.reader(f, delimiter="\t"); samples = next(r)[1:]
    for row in r:
        gc[row[0]] = np.array([float(x) for x in row[1:]])
genera = list(gc); gi = {g: i for i, g in enumerate(genera)}
G = np.vstack([gc[k] for k in genera])

IDX = ['Sebum', 'Moisture', 'Age', 'Burnish', 'Desquamation', 'Temperature', 'Redness']
X = {k: np.array([float(md[s][k]) for s in samples]) for k in IDX}
grp = np.array([md[s]['Group'] for s in samples])

W("=" * 78)
W("14  임상 지표 7개 전수 검정 — 정말 두 축뿐인가")
W(f"    시료 {len(samples)}개 (NS {int((grp=='NS').sum())} / SS {int((grp=='SS').sum())}) · 순열 {NPERM}회 · 난수 {SEED}")
W("    새로 검정하는 것: Burnish(윤기) · Desquamation(각질) · Temperature · Redness")
W("=" * 78)


def shannon(v):
    p = v[v > 0] / v.sum(); return float(-(p * np.log(p)).sum())


def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d: return v.astype(float)
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False),
                       minlength=len(v)).astype(float)


depth = int(G.sum(axis=0).min())
R = np.column_stack([rarefy(G[:, j], depth) for j in range(G.shape[1])])
rel = R / R.sum(axis=0)
sh = np.array([shannon(R[:, j]) for j in range(R.shape[1])])
dom = rel.max(axis=0)
obs_n = (R > 0).sum(axis=0).astype(float)
lr_SL = np.log((R[gi['Staphylococcus']] + .5) / (R[gi['Lawsonella']] + .5))
lr_CL = np.log((R[gi['Cutibacterium']] + .5) / (R[gi['Lawsonella']] + .5))
W(f"\n  공통 깊이 {depth:,} reads")


def pcor(x, y, z=None, n=NPERM):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    if z is None:
        a, b = rx - rx.mean(), ry - ry.mean()
    else:
        Z = np.c_[np.ones(len(x)), stats.rankdata(z)]
        f = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
        a, b = f(rx), f(ry)
    o = np.corrcoef(a, b)[0, 1]
    null = np.array([np.corrcoef(rng.permutation(a), b)[0, 1] for _ in range(n)])
    return o, (np.sum(np.abs(null) >= abs(o)) + 1) / (n + 1)


def bh(ps):
    ps = np.asarray(ps, float); m = len(ps); o = np.argsort(ps)
    q = np.empty(m); prev = 1.0
    for rank, i in enumerate(o[::-1]):
        k = m - rank
        prev = min(prev, ps[i] * m / k); q[i] = prev
    return q


# ------------------------------------------------------------------ A. 지표 7개 x 결과 5종
TARGETS = [("Shannon", sh), ("Observed", obs_n), ("우점도", dom),
           ("log(Staph/Laws)", lr_SL), ("log(Cuti/Laws)", lr_CL)]
W("\n### A. 지표 7개 × 군집 지표 5종 (나이 통제, 지표군 안에서 BH-FDR)")
W("    ※ Age 행만 통제 없이(자기 자신이므로)")
res = {}
for tname, y in TARGETS:
    ps, rs = [], []
    for k in IDX:
        r, p = pcor(X[k], y, None if k == 'Age' else X['Age'])
        rs.append(r); ps.append(p)
    qs = bh(ps)
    res[tname] = dict(zip(IDX, zip(rs, ps, qs)))
    W(f"\n  ── {tname}")
    W(f"     {'지표':<14}{'rho':>8}{'p':>9}{'q(BH)':>9}   판정")
    for k, r, p, q in sorted(zip(IDX, rs, ps, qs), key=lambda t: t[2]):
        mark = "**유의**" if q < 0.05 else ("경향" if q < 0.10 else "")
        W(f"     {k:<14}{r:>+8.3f}{p:>9.4f}{q:>9.4f}   {mark}")

# ------------------------------------------------------------------ B. 새 지표 4개 요약
W("\n### B. 새로 검정한 4개는 무엇을 설명하나")
NEW = ['Burnish', 'Desquamation', 'Temperature', 'Redness']
W(f"  {'지표':<14}" + "".join(f"{t[:13]:>15}" for t, _ in TARGETS))
for k in NEW + ['Sebum', 'Moisture']:
    row = f"  {k:<14}"
    for tname, _ in TARGETS:
        r, p, q = res[tname][k]
        row += f"{f'{r:+.2f}' + ('*' if q < 0.05 else ''):>15}"
    W(row)
W("\n  * = BH-FDR q < 0.05")
hits = [(t, k) for t in res for k in NEW if res[t][k][2] < 0.05]
W(f"\n  새 지표 4개가 q<0.05 로 잡힌 건수: {len(hits)}건 {hits if hits else ''}")

# ------------------------------------------------------------------ C. 새 지표가 피지·수분과 겹치나
W("\n### C. 새 지표 4개가 피지·수분과 얼마나 겹치나 (나이 통제 편상관)")
W(f"  {'지표':<14}{'~ Sebum':>12}{'~ Moisture':>12}")
for k in NEW:
    a, pa = pcor(X[k], X['Sebum'], X['Age'])
    b, pb = pcor(X[k], X['Moisture'], X['Age'])
    W(f"  {k:<14}{f'{a:+.3f}' + ('*' if pa < .05 else ''):>12}{f'{b:+.3f}' + ('*' if pb < .05 else ''):>12}")
W("  * = p<0.05 (보정 전)")

# ------------------------------------------------------------------ D. 군 차이
W("\n### D. 지표별 NS vs SS 군 차이 (Mann-Whitney, BH-FDR)")
ps, us = [], []
for k in IDX:
    u = stats.mannwhitneyu(X[k][grp == 'NS'], X[k][grp == 'SS'])
    us.append(u.statistic); ps.append(u.pvalue)
qs = bh(ps)
W(f"  {'지표':<14}{'NS 중앙':>10}{'SS 중앙':>10}{'p':>10}{'q(BH)':>9}")
for k, p, q in sorted(zip(IDX, ps, qs), key=lambda t: t[1]):
    W(f"  {k:<14}{np.median(X[k][grp=='NS']):>10.2f}{np.median(X[k][grp=='SS']):>10.2f}{p:>10.4f}{q:>9.4f}"
      + ("   **유의**" if q < 0.05 else ""))

# ------------------------------------------------------------------ E. 조성 전체 설명력
W("\n### E. 조성 전체(Bray-Curtis)에 대한 각 지표의 단독 설명력 — PERMANOVA 유사 R²")
from scipy.spatial.distance import braycurtis
n = len(samples)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        D[i, j] = D[j, i] = braycurtis(rel[:, i], rel[:, j])
A = -0.5 * D ** 2
Cn = np.eye(n) - np.ones((n, n)) / n
Gm = Cn @ A @ Cn
tot = np.trace(Gm)


def r2_of(v):
    v = stats.rankdata(v).astype(float)
    v = (v - v.mean()) / v.std()
    H = np.outer(v, v) / (v @ v)
    return np.trace(H @ Gm @ H) / tot if tot > 0 else np.nan


ps, r2s = [], []
for k in IDX:
    o = r2_of(X[k])
    null = np.array([r2_of(rng.permutation(X[k])) for _ in range(999)])
    r2s.append(o); ps.append((np.sum(null >= o) + 1) / 1000)
qs = bh(ps)
W(f"  {'지표':<14}{'R²':>9}{'p':>9}{'q(BH)':>9}")
for k, r2, p, q in sorted(zip(IDX, r2s, ps, qs), key=lambda t: -t[1]):
    W(f"  {k:<14}{r2:>9.4f}{p:>9.4f}{q:>9.4f}" + ("   **유의**" if q < 0.05 else ""))

W("\n" + "=" * 78)
W("### 판정")
if not hits:
    W("  새로 검정한 4개(Burnish·Desquamation·Temperature·Redness) 중 어느 것도")
    W("  군집 지표 5종에 대해 q<0.05 로 잡히지 않았다.")
    W("  -> '임상 지표 일곱 개를 전부 검정했고 피지와 수분 둘만 군집을 구조화한다' 로 쓸 수 있다.")
    W("     이는 두 지표만 골라 본 것보다 강한 주장이다.")
else:
    W(f"  ⚠️ 새 지표에서 {len(hits)}건이 잡혔다. 두 축 프레임을 다시 봐야 한다.")

open(f"{ROOT}/output/14_all_seven.txt", "w").write("\n".join(rep))
print(f"\n저장: {ROOT}/output/14_all_seven.txt")
