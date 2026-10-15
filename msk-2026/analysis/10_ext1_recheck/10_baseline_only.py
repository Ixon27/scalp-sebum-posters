#!/usr/bin/env python
"""
10  외부 코호트 PRJNA1465550 재검 — baseline / treatment 를 갈라서 본다.

배경(2026-08-26 발견):
  이 코호트 시료명 SD##.LS|NLS.1|2 의 접미 .1/.2 를 어제(08-25) 는 '기술 반복'으로 읽고
  사람×부위 단위로 **합산**했다. 그러나 ENA 의 experiment_title 을 보면
      .1 = baseline,  .2 = treatment
  이고 242 run 전부에서 예외 0건으로 일치한다. 피험자는 A군 18명 / K군 22명으로 고정 배정.
  즉 어제 결과는 처치 전후 시료를 섞은 위에 서 있다.

이 스크립트가 하는 일
  (A) baseline 만 / treatment 만 / 어제처럼 합산 — 세 갈래로 07 의 두 검정을 다시 돌린다.
  (B) 처치가 실제로 미생물을 바꿨는지 짝 검정한다. 안 바꿨으면 어제 결론은 살아남는다.

결론을 미리 정하지 않는다. 세 갈래가 어긋나면 어긋난 대로 적는다.
"""
import re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
SEED, NPERM = 20260825, 9999
rng = np.random.default_rng(SEED)

ROOT = Path(__file__).resolve().parents[2]   # = msk-2026/
BASE = ROOT / "data" / "PRJNA1465550"   # QIIME2 export_table/ · export_taxonomy/
OUT = ROOT / "output" / "ext1_baseline"
OUT.mkdir(parents=True, exist_ok=True)
HERE = Path(__file__).resolve().parent   # runs_group.tsv 위치

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


def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d:
        return v.astype(float)
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False),
                       minlength=len(v)).astype(float)


# ------------------------------------------------------------------ 자료
tsv = BASE / "export_table" / "feature-table.tsv"
tax_f = BASE / "export_taxonomy" / "taxonomy.tsv"
if not tsv.exists():
    sys.exit(f"없음: {tsv}")
counts = pd.read_csv(tsv, sep="\t", skiprows=1, index_col=0)
tax = pd.read_csv(tax_f, sep="\t", index_col=0)
G = to_genus(counts, tax)

meta = [(s, m.group(1), m.group(2), int(m.group(3))) for s in G.columns
        if (m := re.match(r'^(SD\d+)\.(LS|NLS)\.(\d+)$', s))]
E = pd.DataFrame(meta, columns=['sample', 'subject', 'site', 'rep']).set_index('sample')
G = G[E.index]

W("=" * 78)
W("10  외부 코호트 재검 — .1/.2 는 기술 반복이 아니라 baseline/treatment 였다")
W(f"    시료 {G.shape[1]}개 · 피험자 {E.subject.nunique()}명 · 속 {G.shape[0]}개 · 순열 {NPERM}회 · seed {SEED}")
W("=" * 78)
W(f"\n  rep 분포: {dict(E.rep.value_counts().sort_index())}   (.1=baseline, .2=treatment)")
W(f"  site 분포: {dict(E.site.value_counts())}")


def metrics(sub_G, P):
    depth = int(sub_G.sum(axis=0).min())
    pr = pd.DataFrame({c: rarefy(sub_G[c].values, depth) for c in sub_G.columns}, index=sub_G.index)
    M = pd.DataFrame({
        'lipid_SL': -np.log((pr.loc['Staphylococcus'] + .5) / (pr.loc['Lawsonella'] + .5)),
        'lipid_CL': -np.log((pr.loc['Cutibacterium'] + .5) / (pr.loc['Lawsonella'] + .5)),
        'shannon': [shannon(pr[c].values) for c in pr.columns],
        'observed': (pr > 0).sum(axis=0).astype(float),
        'dominance': pr.max(axis=0) / pr.sum(axis=0) * 100,
    }, index=pr.columns).join(P)
    return M, depth


def icc_block(Mp, tag):
    W(f"\n  [{tag}] 급내상관 — 사람 간 분산인가 부위 간 분산인가")
    W(f"    {'지표':<24}{'ICC':>8}{'사람간 SD':>11}{'부위간 SD':>11}{'p(사람효과)':>12}")
    icc = {}
    for k in ['lipid_SL', 'lipid_CL', 'shannon', 'observed', 'dominance']:
        a = Mp[Mp.site == 'NLS'].set_index('subject')[k].sort_index()
        b = Mp[Mp.site == 'LS'].set_index('subject')[k].sort_index()
        a, b = a.align(b, join='inner')
        subj_mean = (a + b) / 2
        var_b = np.var(subj_mean, ddof=1)
        var_w = np.mean((b - a) ** 2) / 2
        v = max(var_b - var_w / 2, 0.0)
        ic = v / (v + var_w) if (v + var_w) > 0 else np.nan
        vals = np.concatenate([a.values, b.values])
        null = [np.var(rng.permutation(vals).reshape(2, -1).mean(axis=0), ddof=1) for _ in range(999)]
        pv = (np.sum(np.array(null) >= var_b) + 1) / 1000
        icc[k] = ic
        W(f"    {k:<24}{ic:>8.3f}{np.sqrt(v):>11.3f}{np.sqrt(var_w):>11.3f}{pv:>12.4f}")
    hi = np.mean([icc['lipid_SL'], icc['lipid_CL']])
    lo = np.mean([icc['shannon'], icc['observed']])
    W(f"    지질 축 평균 {hi:.3f} vs 다양성 평균 {lo:.3f}  ->  {'가설과 일치' if hi > lo else '가설과 불일치'}")
    return hi, lo


def ortho_block(Mp, tag):
    W(f"\n  [{tag}] 두 축 직교성 — 피험자 평균으로")
    sub = Mp.groupby('subject')[['lipid_SL', 'lipid_CL', 'shannon', 'observed']].mean()
    res = {}
    for k in ['lipid_SL', 'lipid_CL']:
        for d in ['shannon', 'observed']:
            r, p = perm_rho(sub[k].values, sub[d].values)
            res[f"{k}~{d}"] = (r, p)
            flag = "유의" if p < 0.05 else "비유의"
            W(f"    {k} ~ {d:<9} rho = {r:+.3f}  p = {p:.4f}   {flag}   (n={len(sub)})")
    return res


# ------------------------------------------------------------------ (A) 세 갈래
strata = {
    "baseline만 (.1)":  E[E.rep == 1].index,
    "treatment만 (.2)": E[E.rep == 2].index,
    "합산 (어제 방식)":   None,
}
summary = {}
for tag, idx in strata.items():
    W("\n" + "-" * 78)
    W(f"### {tag}")
    if idx is None:
        key = E.subject + "|" + E.site
        sub_G = G.T.groupby(key.values).sum().T
        P = pd.DataFrame([c.split("|") for c in sub_G.columns],
                         columns=['subject', 'site'], index=sub_G.columns)
    else:
        sub_G = G[idx]
        P = E.loc[idx, ['subject', 'site']]
    paired = sorted(s for s in P.subject.unique() if {'LS', 'NLS'} <= set(P.site[P.subject == s]))
    M, depth = metrics(sub_G, P)
    Mp = M[M.subject.isin(paired)]
    W(f"  단위 {sub_G.shape[1]}개 · 짝 피험자 {len(paired)}명 · 공통 깊이 {depth:,} reads")
    hi, lo = icc_block(Mp, tag)
    res = ortho_block(Mp, tag)
    summary[tag] = dict(n=len(paired), depth=depth, icc_hi=hi, icc_lo=lo, ortho=res)

# ------------------------------------------------------------------ (B) 처치 효과
W("\n" + "=" * 78)
W("### (B) 처치가 실제로 미생물을 바꿨나 — 같은 사람·같은 부위의 baseline vs treatment 짝 검정")
W("    안 바꿨다면 어제의 합산은 결과적으로 무해하다. 바꿨다면 어제 수치는 흔들린다.")
key_bt = E.subject + "|" + E.site
depth_all = int(G.sum(axis=0).min())
pr_all = pd.DataFrame({c: rarefy(G[c].values, depth_all) for c in G.columns}, index=G.index)
Mall = pd.DataFrame({
    'lipid_SL': -np.log((pr_all.loc['Staphylococcus'] + .5) / (pr_all.loc['Lawsonella'] + .5)),
    'lipid_CL': -np.log((pr_all.loc['Cutibacterium'] + .5) / (pr_all.loc['Lawsonella'] + .5)),
    'shannon': [shannon(pr_all[c].values) for c in pr_all.columns],
    'observed': (pr_all > 0).sum(axis=0).astype(float),
    'dominance': pr_all.max(axis=0) / pr_all.sum(axis=0) * 100,
}, index=pr_all.columns).join(E)
W(f"    전 시료 공통 깊이 {depth_all:,} reads")
W(f"\n    {'지표':<24}{'baseline':>11}{'treatment':>11}{'차이':>10}{'p(Wilcoxon)':>13}{'짝 수':>7}")
for k in ['lipid_SL', 'lipid_CL', 'shannon', 'observed', 'dominance']:
    a = Mall[Mall.rep == 1].set_index(key_bt[Mall[Mall.rep == 1].index])[k]
    b = Mall[Mall.rep == 2].set_index(key_bt[Mall[Mall.rep == 2].index])[k]
    a, b = a.align(b, join='inner')
    if len(a) < 5:
        W(f"    {k:<24} 짝 부족 (n={len(a)})"); continue
    st = stats.wilcoxon(a.values, b.values)
    W(f"    {k:<24}{np.median(a):>11.3f}{np.median(b):>11.3f}{np.median(b-a):>10.3f}{st.pvalue:>13.4f}{len(a):>7}")

W("\n### (C) 군(A/K) 별로도 갈리나 — baseline 만, 군 안에서 직교성")
grp_map = pd.read_csv(HERE / "runs_group.tsv", sep="\t", index_col=0)['grp'] \
    if (HERE / "runs_group.tsv").exists() else None
if grp_map is None:
    W("    runs_group.tsv 없음 — 군별 분해 생략")
else:
    Eb = E[E.rep == 1].copy()
    Eb['grp'] = Eb.subject.map(grp_map)
    for g in sorted(Eb.grp.dropna().unique()):
        idx = Eb[Eb.grp == g].index
        P = Eb.loc[idx, ['subject', 'site']]
        paired = sorted(s for s in P.subject.unique() if {'LS', 'NLS'} <= set(P.site[P.subject == s]))
        if len(paired) < 6:
            W(f"    {g}군: 짝 피험자 {len(paired)}명 — 너무 적어 생략"); continue
        M, depth = metrics(G[idx], P)
        ortho_block(M[M.subject.isin(paired)], f"{g}군 baseline")

# ------------------------------------------------------------------ 요약
W("\n" + "=" * 78)
W("### 요약 — 어제 결론이 살아남는가")
W(f"    {'갈래':<20}{'짝 n':>6}{'지질ICC':>9}{'다양성ICC':>10}   직교성(유의 개수/4)")
for tag, s in summary.items():
    nsig = sum(1 for r, p in s['ortho'].values() if p < 0.05)
    W(f"    {tag:<20}{s['n']:>6}{s['icc_hi']:>9.3f}{s['icc_lo']:>10.3f}   {nsig}/4")
W("\n    직교성이 재현됐다면 '유의 개수'가 0/4 이어야 한다.")

(OUT / "10_baseline_only.txt").write_text("\n".join(rep), encoding="utf-8")
print(f"\n저장: {OUT/'10_baseline_only.txt'}")
