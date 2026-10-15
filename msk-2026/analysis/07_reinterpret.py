#!/usr/bin/env python
"""
07  외부 검증 실패를 해석한다 — 사후 해석이므로 검정을 붙인다.

배경: 우리 코호트에서 정의한 지질 축(로그비)이 외부 코호트의 병변/비병변 대비에서
      재현되지 않았고, 음성대조로 걸어둔 다양성 쪽이 오히려 유의하게 움직였다.

가능한 설명 두 가지를 구분한다.
  (H1) 우리 축 분리가 틀렸다.
  (H2) 검정이 어긋났다 — 지질 축은 '사람 단위'(피지는 개인 형질) 성질인데,
       같은 사람의 병변/비병변은 '부위 단위'(국소 염증) 대비라 지질 축이 움직일 이유가 없다.

H2가 맞다면 두 가지가 성립해야 한다. 둘 다 검정한다.
  (검정 1) 급내상관 — 지질 축은 사람 간 분산이 지배하고, 다양성은 부위 간 분산이 크다.
  (검정 2) 축 독립성이 외부 코호트에서도 재현된다 — 임상 자료 없이 미생물만으로 확인 가능.
           우리 코호트에서 지질 균형과 다양성은 직교했다. 외부에서도 직교하면 두 축 구조가 재현된 것이다.

사후 해석임을 잊지 않는다. 검정을 통과해도 '가설과 모순되지 않는다'까지가 정직한 문장이다.
"""
import os, re, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
SEED, NPERM = 20260825, 9999
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
rng = np.random.default_rng(SEED)
rep = []


def W(*a):
    s = "".join(str(x) for x in a); rep.append(s); print(s, flush=True)


def rank_genus(t):
    for p in (x.strip() for x in t.split(';')):
        if p.startswith('g__'):
            return p[3:]
    return ''


def to_genus(counts, tax):
    gg = tax['Taxon'].map(rank_genus); gg = gg[gg != '']
    idx = counts.index.intersection(gg.index)
    return counts.loc[idx].join(gg.rename('genus')).groupby('genus').sum()


def shannon(v):
    p = v[v > 0] / v.sum()
    return float(-(p * np.log(p)).sum())


def perm_rho(x, y, n=NPERM):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    obs = np.corrcoef(rx, ry)[0, 1]
    null = np.array([np.corrcoef(rng.permutation(rx), ry)[0, 1] for _ in range(n)])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (n + 1)


# ------------------------------------------------------------------ 자료
tax_own = pd.read_csv(ROOT / "data/external/own_taxonomy_silva138.tsv", sep="\t", index_col=0)
asv_own = pd.read_csv(ROOT / "data/own_asv_table.tsv", sep="\t", index_col=0)
md = pd.read_csv(ROOT / "output/01_metadata.tsv", sep="\t", index_col=0)
own = to_genus(asv_own, tax_own)[md.index]

ext_counts = pd.read_csv(ROOT / "data/external/feature-table.tsv", sep="\t", skiprows=1, index_col=0)
ext_tax = pd.read_csv(ROOT / "data/external/taxonomy.tsv", sep="\t", index_col=0)
ext = to_genus(ext_counts, ext_tax)

meta = [(s, m.group(1), m.group(2)) for s in ext.columns
        if (m := re.match(r'^(SD\d+)\.(LS|NLS)\.\d+$', s))]
E = pd.DataFrame(meta, columns=['sample', 'subject', 'site']).set_index('sample')
ext = ext[E.index]
pool = ext.T.groupby((E.subject + "|" + E.site).values).sum().T
P = pd.DataFrame([c.split("|") for c in pool.columns], columns=['subject', 'site'], index=pool.columns)

depth = int(pool.sum(axis=0).min())


def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d:
        return v.astype(float)
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False),
                       minlength=len(v)).astype(float)


pr = pd.DataFrame({c: rarefy(pool[c].values, depth) for c in pool.columns}, index=pool.index)

M = pd.DataFrame({
    'lipid_SL': -np.log((pr.loc['Staphylococcus'] + .5) / (pr.loc['Lawsonella'] + .5)),
    'lipid_CL': -np.log((pr.loc['Cutibacterium'] + .5) / (pr.loc['Lawsonella'] + .5)),
    'shannon': [shannon(pr[c].values) for c in pr.columns],
    'observed': (pr > 0).sum(axis=0).astype(float),
    'dominance': pr.max(axis=0) / pr.sum(axis=0) * 100,
}, index=pr.columns).join(P)
paired = sorted(s for s in P.subject.unique() if {'LS', 'NLS'} <= set(P.site[P.subject == s]))
Mp = M[M.subject.isin(paired)]

W("=" * 78)
W("07  외부 검증 실패의 해석 — 사후 가설을 검정한다")
W(f"    외부 코호트 짝 피험자 {len(paired)}명, 공통 깊이 {depth:,} reads, 순열 {NPERM}회")
W("=" * 78)

# ============================================================= 검정 1
W("\n### 검정 1 — 급내상관(ICC): 사람 간 분산인가, 부위 간 분산인가")
W("  ICC = 사람 간 분산 / 전체 분산. 1에 가까울수록 '사람의 형질', 0에 가까울수록 '부위의 상태'.")
W(f"\n  {'지표':<26}{'ICC':>8}{'사람간 SD':>11}{'부위간 SD':>11}{'p(사람 효과)':>13}")
icc = {}
for k in ['lipid_SL', 'lipid_CL', 'shannon', 'observed', 'dominance']:
    a = Mp[Mp.site == 'NLS'].set_index('subject')[k].sort_index()
    b = Mp[Mp.site == 'LS'].set_index('subject')[k].sort_index()
    a, b = a.align(b, join='inner')
    subj_mean = (a + b) / 2
    var_b = np.var(subj_mean, ddof=1)                      # 사람 간
    var_w = np.mean((b - a) ** 2) / 2                      # 부위 간 (짝 차이로부터)
    v = var_b - var_w / 2
    v = max(v, 0.0)
    ic = v / (v + var_w) if (v + var_w) > 0 else np.nan
    # 사람 효과 유의성: 사람 라벨을 섞어 사람 간 분산 귀무분포를 만든다
    vals = np.concatenate([a.values, b.values])
    null = []
    for _ in range(999):
        p2 = rng.permutation(vals).reshape(2, -1)
        null.append(np.var(p2.mean(axis=0), ddof=1))
    pv = (np.sum(np.array(null) >= var_b) + 1) / 1000
    icc[k] = ic
    W(f"  {k:<26}{ic:>8.3f}{np.sqrt(v):>11.3f}{np.sqrt(var_w):>11.3f}{pv:>13.4f}")
W("\n  가설: 지질 축(lipid_*)의 ICC 가 다양성(shannon·observed)보다 높아야 한다.")
hi = np.mean([icc['lipid_SL'], icc['lipid_CL']])
lo = np.mean([icc['shannon'], icc['observed']])
W(f"  지질 축 평균 ICC {hi:.3f}  vs  다양성 평균 ICC {lo:.3f}  ->  "
  f"{'가설과 일치' if hi > lo else '가설과 불일치'}")

# ============================================================= 검정 2
W("\n### 검정 2 — 두 축의 독립성이 외부 코호트에서도 재현되나")
W("  임상 자료가 없어도 확인할 수 있다. 우리 코호트에서 지질 균형과 다양성은 직교했다.")

W("\n  [우리 코호트, n=27] 나이 통제 안 한 단순 순위상관")
own_r = own / own.sum(axis=0)
own_sh = np.array([shannon(own[c].values) for c in own.columns])
own_lr = -np.log((own.loc['Staphylococcus'] + .5) / (own.loc['Lawsonella'] + .5)).values
r0, p0 = perm_rho(own_lr, own_sh)
W(f"    지질 균형 ~ Shannon    rho = {r0:+.3f}  p = {p0:.4f}")

W(f"\n  [외부 코호트, 피험자 {len(paired)}명 평균] 같은 상관")
sub = Mp.groupby('subject')[['lipid_SL', 'lipid_CL', 'shannon', 'observed', 'dominance']].mean()
for k in ['lipid_SL', 'lipid_CL']:
    for d in ['shannon', 'observed']:
        r, p = perm_rho(sub[k].values, sub[d].values)
        W(f"    {k} ~ {d:<10} rho = {r:+.3f}  p = {p:.4f}")
W("\n  가설: 외부에서도 지질 균형과 다양성이 서로 독립(비유의)이어야 두 축 구조가 재현된 것이다.")

# ============================================================= 코호트 차이
W("\n### 참고 — 두 코호트의 군집 구조가 얼마나 다른가")
W(f"  {'속':<20}{'우리 평균 RA':>14}{'외부 평균 RA':>14}")
ext_ra = (pr / pr.sum(axis=0) * 100).mean(axis=1)
own_ra = (own / own.sum(axis=0) * 100).mean(axis=1)
for g in ['Staphylococcus', 'Cutibacterium', 'Lawsonella', 'Streptococcus', 'Corynebacterium']:
    W(f"  {g:<20}{own_ra.get(g,0):>13.2f}%{ext_ra.get(g,0):>13.2f}%")
W(f"  {'Shannon 중앙':<20}{np.median(own_sh):>13.2f} {np.median(M.shannon):>13.2f}")
W("\n  두 코호트는 우점 속 순서 자체가 뒤집혀 있다. 절대값 비교는 하지 않고 방향과 짝 차이만 본다.")

(OUT / "07_reinterpret.txt").write_text("\n".join(rep), encoding="utf-8")
print("\n저장: output/07_reinterpret.txt")
