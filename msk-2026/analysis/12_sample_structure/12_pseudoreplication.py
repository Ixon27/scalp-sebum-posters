#!/usr/bin/env python
"""
12  우리 27개 시료가 정말 27명인가 — 유사반복(pseudoreplication) 검정.

발단 (2026-08-26): 나이가 27개 시료에 13종밖에 없고, 같은 나이끼리 임상값이 지나치게 닮았다.
  NS01/NS02 (둘 다 24세) : Moisture 2.76 = 2.76,  Sebum 21 = 21,  Redness 106.33 = 106.33  (완전 일치)
  NS06/07/08 (58세)      : Moisture 6.44/6.40/7.60,  Temperature 32.2/32.2/32.8
  NS09/10/11 (61세)      : Redness 175.00/175.67/167.00

가설
  H0 : 27명이 서로 다른 사람이고 나이 일치는 우연이다.
  H1 : 실제 사람은 13명이고 한 사람에게서 2~3개 시료를 얻었다 (유사반복).

H1 이 맞다면 **미생물 조성도** 같은 나이 묶음 안에서 훨씬 닮아야 한다.
판정 기준선은 오늘 외부 코호트에서 실측한 값을 쓴다 (11_dup_runs):
  같은 사람 다른 시기 Bray-Curtis 중앙 0.354  <  다른 사람 중앙 0.650

이게 사실이면 유효 표본이 27 -> 13 으로 줄고 모든 p 값이 부풀려진 것이다. 반드시 확인한다.
"""
import csv, itertools, collections, sys
import numpy as np
from scipy import stats
from scipy.spatial.distance import braycurtis

SEED, NPERM = 20260825, 9999
rng = np.random.default_rng(SEED)
ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
rep = []
def W(*a):
    s = "".join(str(x) for x in a); rep.append(s); print(s, flush=True)


# ---------------------------------------------------------------- 자료
md = {}
with open(f"{ROOT}/output/01_metadata.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    idc = r.fieldnames[0]
    for row in r:
        md[row[idc]] = row

gc = {}
with open(f"{ROOT}/output/01_genus_counts.tsv") as f:
    r = csv.reader(f, delimiter="\t")
    hdr = next(r)
    samples = hdr[1:]
    for row in r:
        gc[row[0]] = np.array([float(x) for x in row[1:]])
G = np.vstack([gc[k] for k in gc])          # 속 x 시료
genera = list(gc)
W("=" * 78)
W("12  27개 시료가 27명인가 — 유사반복 검정")
W("=" * 78)
W(f"\n  속 {G.shape[0]}개 · 시료 {G.shape[1]}개 · 시료명 {samples[:3]} …")

ages = {s: float(md[s]['Age']) for s in samples}
by_age = collections.defaultdict(list)
for s in samples:
    by_age[ages[s]].append(s)
W(f"  나이 고유값 {len(by_age)}종 / 시료 {len(samples)}개")
W(f"  나이별 시료 수: {dict(collections.Counter(len(v) for v in by_age.values()))}  (2 이상이면 묶음)")
W("\n  나이 묶음:")
for a in sorted(by_age):
    v = by_age[a]
    W(f"    {a:>5.0f}세 : {v}")

# ---------------------------------------------------------------- 임상값 일치
W("\n### A. 임상값이 나이 묶음 안에서 얼마나 겹치나")
clin = ['Moisture', 'Sebum', 'Temperature', 'Redness', 'Burnish', 'Desquamation']
W(f"  {'변수':<14}{'묶음 내 완전일치 쌍':>18}{'전체 쌍 중 일치':>16}")
within_pairs = [p for a in by_age for p in itertools.combinations(by_age[a], 2)]
all_pairs = list(itertools.combinations(samples, 2))
for c in clin:
    w = sum(1 for x, y in within_pairs if md[x][c] == md[y][c])
    t = sum(1 for x, y in all_pairs if md[x][c] == md[y][c])
    W(f"  {c:<14}{f'{w}/{len(within_pairs)}':>18}{f'{t}/{len(all_pairs)}':>16}")
W(f"\n  묶음 내 쌍 {len(within_pairs)}개, 전체 쌍 {len(all_pairs)}개")

# ---------------------------------------------------------------- 미생물 거리
W("\n### B. 미생물 조성 — 나이 묶음 안이 밖보다 닮았나 (핵심 검정)")
rel = G / G.sum(axis=0, keepdims=True)
idx = {s: i for i, s in enumerate(samples)}
D = {}
for x, y in all_pairs:
    D[(x, y)] = braycurtis(rel[:, idx[x]], rel[:, idx[y]])

win = np.array([D[p] if p in D else D[(p[1], p[0])] for p in within_pairs])
wset = set(within_pairs) | {(b, a) for a, b in within_pairs}
btw = np.array([d for p, d in D.items() if p not in wset])

W(f"\n  {'비교':<26}{'n':>6}{'중앙 BC':>10}{'25%':>8}{'75%':>8}")
W(f"  {'같은 나이 묶음 안':<26}{len(win):>6}{np.median(win):>10.3f}{np.percentile(win,25):>8.3f}{np.percentile(win,75):>8.3f}")
W(f"  {'다른 나이':<26}{len(btw):>6}{np.median(btw):>10.3f}{np.percentile(btw,25):>8.3f}{np.percentile(btw,75):>8.3f}")

obs = np.median(btw) - np.median(win)
lab = np.array([ages[s] for s in samples])
null = []
for _ in range(NPERM):
    perm = rng.permutation(lab)
    pa = {s: perm[i] for i, s in enumerate(samples)}
    w2, b2 = [], []
    for (x, y), d in D.items():
        (w2 if pa[x] == pa[y] else b2).append(d)
    null.append(np.median(b2) - np.median(w2) if w2 else 0.0)
null = np.array(null)
p = (np.sum(null >= obs) + 1) / (NPERM + 1)
W(f"\n  중앙값 차이 {obs:+.3f}  순열 p = {p:.4f}   (나이 라벨을 섞어 만든 귀무분포, {NPERM}회)")

W("\n  [기준선] 오늘 외부 코호트에서 실측한 값 (11_dup_runs)")
W("    같은 사람·다른 시기  BC 중앙 0.354")
W("    같은 사람·다른 부위  BC 중앙 0.403")
W("    다른 사람            BC 중앙 0.650")
W(f"\n  우리 자료의 '같은 나이 묶음 안' = {np.median(win):.3f}")
if np.median(win) <= 0.45:
    W("    -> 외부 코호트의 '같은 사람' 수준. 유사반복 의심이 강해진다.")
elif np.median(win) >= 0.60:
    W("    -> 외부 코호트의 '다른 사람' 수준. 서로 다른 사람이라는 쪽과 맞는다.")
else:
    W("    -> 중간. 이것만으로 판정하지 않는다.")

# ---------------------------------------------------------------- 묶음별 상세
W("\n### C. 묶음별 상세 — 어느 묶음이 특히 닮았나")
W(f"  {'나이':<7}{'n':>3}{'묶음내 중앙 BC':>15}{'임상 완전일치 변수':>22}")
for a in sorted(by_age):
    v = by_age[a]
    if len(v) < 2:
        continue
    ds = [D[p] if p in D else D[(p[1], p[0])] for p in itertools.combinations(v, 2)]
    same = [c for c in clin if len({md[s][c] for s in v}) == 1]
    W(f"  {a:<7.0f}{len(v):>3}{np.median(ds):>15.3f}{(','.join(same) if same else '-'):>22}")

# ---------------------------------------------------------------- 군 구조 확인
W("\n### D. 나이 묶음이 군(NS/SS) 을 넘나드나")
mixed = [a for a in by_age if len({md[s]['Group'] for s in by_age[a]}) > 1]
W(f"  NS 와 SS 가 섞인 나이 묶음: {mixed if mixed else '없음'}")
W("  (섞이지 않는다면 '같은 사람' 가설과 모순되지 않는다 — 한 사람은 한 군이므로)")

open(f"{ROOT}/output/12_pseudoreplication.txt", "w").write("\n".join(rep))
print(f"\n저장: {ROOT}/output/12_pseudoreplication.txt")
