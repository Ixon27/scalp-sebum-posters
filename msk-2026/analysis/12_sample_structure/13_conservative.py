#!/usr/bin/env python
"""
13  최악의 경우 검정 — 같은 나이·같은 군 시료를 '한 사람' 으로 합치면 결론이 버티는가.

12 에서 나온 것: 같은 나이 묶음 안의 미생물 조성이 지나치게 닮았다(BC 중앙 0.117).
사람 식별자가 없어 유사반복 여부를 자료만으로 확정할 수 없다.

그래서 확정을 기다리지 않고 **하한**을 먼저 잰다.
  가정: 같은 나이 + 같은 군 = 같은 사람 (가장 보수적. 한 사람이 두 군일 수 없으므로 군은 가른다)
  이 가정에서 유효 표본이 27 -> 15 로 준다.

세 갈래를 나란히 낸다.
  (a) 원래대로 시료 27개
  (b) 사람 단위 평균 (나이×군 묶음의 임상값·지표를 평균)
  (c) 사람 단위 합산 (count 를 합쳐 다시 rarefy — 조성 관점의 합침)

결론이 세 갈래 모두에서 같으면 유사반복 여부와 무관하게 안전하다.
어긋나면 사람 식별자 확인이 초록 제출의 선행 조건이 된다.
"""
import csv, collections
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
    r = csv.reader(f, delimiter="\t"); hdr = next(r); samples = hdr[1:]
    for row in r:
        gc[row[0]] = np.array([float(x) for x in row[1:]])
genera = list(gc)
G = np.vstack([gc[k] for k in genera])

KEEP = ['Cutibacterium', 'Lawsonella', 'Staphylococcus', 'Streptococcus']
gi = {g: i for i, g in enumerate(genera)}


def shannon(v):
    p = v[v > 0] / v.sum(); return float(-(p * np.log(p)).sum())


def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d: return v.astype(float)
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False),
                       minlength=len(v)).astype(float)


def pcor_perm(x, y, z, n=NPERM):
    """z 를 통제한 순위 편상관 + 순열검정"""
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    Z = np.c_[np.ones_like(rz, dtype=float), rz]
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    a, b = res(rx), res(ry)
    obs = np.corrcoef(a, b)[0, 1]
    null = np.array([np.corrcoef(rng.permutation(a), b)[0, 1] for _ in range(n)])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (n + 1)


def block(tag, cnt, clin):
    """cnt: 속x단위 count, clin: dict of arrays"""
    depth = int(cnt.sum(axis=0).min())
    R = np.column_stack([rarefy(cnt[:, j], depth) for j in range(cnt.shape[1])])
    rel = R / R.sum(axis=0)
    sh = np.array([shannon(R[:, j]) for j in range(R.shape[1])])
    dom = rel.max(axis=0)
    pc = 0.5
    lr_SL = np.log((R[gi['Staphylococcus']] + pc) / (R[gi['Lawsonella']] + pc))
    lr_CL = np.log((R[gi['Cutibacterium']] + pc) / (R[gi['Lawsonella']] + pc))
    seb, moi, age = clin['Sebum'], clin['Moisture'], clin['Age']
    n = cnt.shape[1]
    W(f"\n  ── {tag}  (단위 {n}개, 공통 깊이 {depth:,} reads)")
    out = {}
    r, p = pcor_perm(seb, moi, age); out['두축독립'] = (r, p)
    W(f"     {'피지 ~ 수분 (나이 통제)':<34}rho={r:+.3f}  p={p:.4f}   {'독립' if p>=0.05 else '독립 아님'}")
    for nm, y in (("수분 ~ Shannon (나이 통제)", sh), ("수분 ~ 우점도 (나이 통제)", dom)):
        r, p = pcor_perm(moi, y, age); out[nm] = (r, p)
        W(f"     {nm:<34}rho={r:+.3f}  p={p:.4f}   {'유의' if p<0.05 else '비유의'}")
    for nm, y in (("피지 ~ Shannon (나이 통제)", sh),):
        r, p = pcor_perm(seb, y, age); out[nm] = (r, p)
        W(f"     {nm:<34}rho={r:+.3f}  p={p:.4f}   {'유의' if p<0.05 else '비유의'}")
    for nm, y in (("피지 ~ log(Staph/Laws) (나이 통제)", lr_SL),
                  ("피지 ~ log(Cuti/Laws) (나이 통제)", lr_CL)):
        r, p = pcor_perm(seb, y, age); out[nm] = (r, p)
        W(f"     {nm:<34}rho={r:+.3f}  p={p:.4f}   {'유의' if p<0.05 else '비유의'}")
    for nm, y in (("수분 ~ log(Staph/Laws) (나이 통제)", lr_SL),):
        r, p = pcor_perm(moi, y, age); out[nm] = (r, p)
        W(f"     {nm:<34}rho={r:+.3f}  p={p:.4f}   {'유의' if p<0.05 else '비유의'}")
    return out


W("=" * 78)
W("13  유사반복을 가정한 최악의 경우 — 결론이 버티는가")
W(f"    순열 {NPERM}회 · 난수 {SEED} · 나이는 전 검정에서 통제")
W("=" * 78)

# (a) 원래
clin_a = {k: np.array([float(md[s][k]) for s in samples]) for k in ('Sebum', 'Moisture', 'Age')}
res_a = block("(a) 시료 단위 — 원래 방식", G, clin_a)

# 묶음 정의: 나이 x 군
key = [(md[s]['Age'], md[s]['Group']) for s in samples]
units = sorted(set(key), key=lambda t: (float(t[0]), t[1]))
W(f"\n\n  묶음 정의: 같은 나이 + 같은 군 = 한 사람  ->  단위 {len(units)}개")
W(f"  {'단위':<14}{'n':>3}  시료")
for u in units:
    v = [s for s, k in zip(samples, key) if k == u]
    W(f"  {u[0]:>3}세 {u[1]:<8}{len(v):>3}  {', '.join(v)}")

# (b) 사람 단위 평균
cols_b, clin_b = [], {k: [] for k in ('Sebum', 'Moisture', 'Age')}
for u in units:
    js = [j for j, k in enumerate(key) if k == u]
    rel_j = G[:, js] / G[:, js].sum(axis=0, keepdims=True)
    cols_b.append(rel_j.mean(axis=1) * 3037)          # 평균 조성을 공통 깊이로 되돌림
    for k in clin_b:
        clin_b[k].append(np.mean([float(md[samples[j]][k]) for j in js]))
Gb = np.column_stack(cols_b)
res_b = block("(b) 사람 단위 — 조성 평균", Gb, {k: np.array(v) for k, v in clin_b.items()})

# (c) 사람 단위 합산
cols_c = []
for u in units:
    js = [j for j, k in enumerate(key) if k == u]
    cols_c.append(G[:, js].sum(axis=1))
Gc = np.column_stack(cols_c)
res_c = block("(c) 사람 단위 — count 합산 후 재rarefy", Gc, {k: np.array(v) for k, v in clin_b.items()})

# ---------------------------------------------------------------- 요약
W("\n" + "=" * 78)
W("### 요약 — 세 갈래 비교")
W(f"  {'검정':<34}{'(a) n=27':>16}{'(b) n=15 평균':>16}{'(c) n=15 합산':>16}")
for k in res_a:
    def f(d):
        r, p = d[k]; return f"{r:+.2f} p={p:.3f}"
    W(f"  {k:<34}{f(res_a):>16}{f(res_b):>16}{f(res_c):>16}")

W("\n  판정 기준: (a) 에서 유의했던 것이 (b)·(c) 에서도 유의하면 유사반복 여부와 무관하게 결론이 선다.")
W("            (a) 에서만 유의하면 그 수치는 유사반복이 만든 것일 수 있다 -> 사람 식별자 확인이 선행 조건.")

open(f"{ROOT}/output/13_conservative.txt", "w").write("\n".join(rep))
print(f"\n저장: {ROOT}/output/13_conservative.txt")
