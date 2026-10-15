#!/usr/bin/env python
"""
19  외부 코호트 PRJNA1465550 — 병변(LS) vs 비병변(NLS) 짝 검정을 baseline(.1) 시료만으로 다시 (2026-10-14).

왜: 포스터 F5A 의 p 6개는 06_server_ml_external.py 의 값인데, 06 은 .1/.2(baseline/treatment)를
    '기술 반복'으로 읽어 합산한 34명 분석이다. 10_baseline_only.py 가 .1/.2 의 정체를 밝힌 뒤
    ICC·직교성은 baseline 30명으로 다시 냈지만, F5A 의 짝 검정만 합산 값으로 남아 있었다 (Codex 검토 지적).
무엇이 같은가: 06 과 같은 속 접기, 공통 깊이 rarefy, 지표 정의, ElasticNet 설정(enet), 공통 속 CLR, 난수 20260825.
무엇이 다른가: 합산하지 않고 rep == 1(baseline) 시료만 쓴다. 비교용으로 합산 판도 같은 코드로 다시 낸다.
실행: python analysis/10_ext1_recheck/19_baseline_paired.py   (msk-2026/ 에서)
"""
import re, inspect, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
SEED = 20260825
ROOT = Path(__file__).resolve().parents[2]   # = msk-2026/
OUT = ROOT / "output"
rep = []


def W(*a):
    s = "".join(str(x) for x in a); rep.append(s); print(s, flush=True)


_ALPHA_KW = 'n_alphas' if 'n_alphas' in inspect.signature(ElasticNetCV).parameters else 'alphas'


def enet():   # 06 과 같은 설정
    return make_pipeline(StandardScaler(),
                         ElasticNetCV(l1_ratio=[.5, .9, 1.0], cv=5, max_iter=2000, tol=1e-3,
                                      selection='random', random_state=SEED, n_jobs=1,
                                      **{_ALPHA_KW: 20}))


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


tax_own = pd.read_csv(ROOT / "data/external/own_taxonomy_silva138.tsv", sep="\t", index_col=0)
asv_own = pd.read_csv(ROOT / "data/own_asv_table.tsv", sep="\t", index_col=0)
md = pd.read_csv(ROOT / "output/01_metadata.tsv", sep="\t", index_col=0)
own = to_genus(asv_own, tax_own)[md.index]
ext_counts = pd.read_csv(ROOT / "data/external/feature-table.tsv", sep="\t", skiprows=1, index_col=0)
ext = to_genus(ext_counts, pd.read_csv(ROOT / "data/external/taxonomy.tsv", sep="\t", index_col=0))

meta = [(s, m.group(1), m.group(2), int(m.group(3))) for s in ext.columns
        if (m := re.match(r'^(SD\d+)\.(LS|NLS)\.(\d+)$', s))]
E = pd.DataFrame(meta, columns=['sample', 'subject', 'site', 'rep']).set_index('sample')
ext = ext[E.index]

LAB = [('lipid_SL', 'log(Laws/Staph)'), ('lipid_CL', 'log(Laws/Cuti)'), ('score_lipid', 'sebum model score'),
       ('shannon', 'Shannon'), ('observed', 'observed genera'), ('dominance', 'dominance'),
       ('score_hydration', 'hydration model score')]


def run(tag, units, P):
    rng = np.random.default_rng(SEED)
    paired = sorted(s for s in P.subject.unique() if {'LS', 'NLS'} <= set(P.site[P.subject == s]))
    depth = int(units.sum(axis=0).min())

    def rarefy(v, d):
        v = v.astype(int)
        if v.sum() <= d:
            return v.astype(float)
        pick = rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False)
        return np.bincount(pick, minlength=len(v)).astype(float)

    R = pd.DataFrame({c: rarefy(units[c].values, depth) for c in units.columns}, index=units.index)
    M = pd.DataFrame({
        'lipid_SL': -np.log((R.loc['Staphylococcus'] + .5) / (R.loc['Lawsonella'] + .5)),
        'lipid_CL': -np.log((R.loc['Cutibacterium'] + .5) / (R.loc['Lawsonella'] + .5)),
        'shannon': [shannon(R[c].values) for c in R.columns],
        'observed': (R > 0).sum(axis=0).astype(float),
        'dominance': R.max(axis=0) / R.sum(axis=0) * 100,
    }, index=R.columns).join(P)
    common = sorted(set(own.index) & set(R.index))

    def clr_of(tab):
        x = np.log(tab.loc[common].astype(float) + 0.5)
        return (x - x.mean(axis=0)).T.values

    Xo, Xe = clr_of(own), clr_of(R)
    for target, col in [('Sebum', 'score_lipid'), ('Moisture', 'score_hydration')]:
        M[col] = enet().fit(Xo, md[target].values.astype(float)).predict(Xe)

    W(f"\n### {tag} — 단위 {units.shape[1]}개, 짝 피험자 {len(paired)}명, 공통 깊이 {depth:,} reads, 공통 속 {len(common)}")
    W(f"  {'지표':<24}{'NLS':>10}{'LS':>10}{'차이 중앙':>11}{'p':>9}{'짝':>5}")
    rows = []
    for k, lab in LAB:
        a = M[(M.site == 'NLS') & M.subject.isin(paired)].set_index('subject')[k].sort_index()
        b = M[(M.site == 'LS') & M.subject.isin(paired)].set_index('subject')[k].sort_index()
        a, b = a.align(b, join='inner')
        p = stats.wilcoxon(b.values, a.values).pvalue
        d = np.median(b.values - a.values)
        W(f"  {lab:<24}{np.median(a):>10.3f}{np.median(b):>10.3f}{d:>+11.3f}{p:>9.4f}{len(a):>5}")
        rows.append(dict(stratum=tag, metric=k, label=lab, n_pairs=len(a), median_NLS=np.median(a),
                         median_LS=np.median(b), median_diff=d, p=p))
    return rows


W("=" * 78)
W("19  PRJNA1465550 병변 vs 비병변 짝 Wilcoxon signed-rank — baseline 만 vs 합산(06 방식)")
W("=" * 78)
rows = []
b = E[E.rep == 1].index
rows += run("baseline (.1)", ext[b], E.loc[b, ['subject', 'site']])
key = E.subject + "|" + E.site
pool = ext.T.groupby(key.values).sum().T
Pp = pd.DataFrame([c.split("|") for c in pool.columns], columns=['subject', 'site'], index=pool.columns)
rows += run("합산 (06 방식, 재현 확인)", pool, Pp)
pd.DataFrame(rows).to_csv(OUT / "19_baseline_paired.tsv", sep="\t", index=False)
(OUT / "19_baseline_paired.txt").write_text("\n".join(rep), encoding="utf-8")
print("저장: output/19_baseline_paired.txt / .tsv")
