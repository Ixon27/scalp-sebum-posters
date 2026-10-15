#!/usr/bin/env python
"""
S2-1(2부)  독립 코호트 외부 검증 — 음성대조가 내장된 예측.

외부 코호트: PRJNA1465550 (지루성 피부염 환자 40명, 같은 사람의 병변 LS / 비병변 NLS 짝,
             16S V3-V4, 우리 데이터와 같은 SILVA 138 분류기로 재분류)

예측 (틀릴 수 있게 만든 가설)
   지질 축 지표는 병변에서 유의하게 높고,   <- 맞아야 함
   수분 축 지표(다양성)는 차이가 없다.       <- 음성대조. 여기까지 같이 움직이면 축 분리가 틀린 것.

짝지은 사람 안에서 비교하므로 개인차·배치 효과가 상당 부분 상쇄된다.
절대값은 코호트 간 비교하지 않는다 — 방향과 짝 차이만 본다.
"""
import warnings, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import ElasticNetCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")

SEED = 20260825
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
EXT = ROOT / "data" / "external"

rep = []; W = rep.append


def collapse_to_genus(counts, tax):
    def rank(s, k):
        for p in (x.strip() for x in s.split(';')):
            if p.startswith(k + '__'):
                return p[3:]
        return ''
    gcol = tax['Taxon'].map(lambda s: rank(s, 'g'))
    gcol = gcol[gcol != '']
    common = counts.index.intersection(gcol.index)
    return counts.loc[common].join(gcol.rename('genus')).groupby('genus').sum()


# ---------------------------------------------------------------- 자료
ext_counts = pd.read_csv(EXT / "feature-table.tsv", sep="\t", skiprows=1, index_col=0)
ext_tax = pd.read_csv(EXT / "taxonomy.tsv", sep="\t", index_col=0)
ext_g = collapse_to_genus(ext_counts, ext_tax)

own_counts = pd.read_csv(OUT / "01_genus_counts.tsv", sep="\t", index_col=0)   # 원 분류기 기준
md = pd.read_csv(OUT / "01_metadata.tsv", sep="\t", index_col=0)
# 우리 데이터를 같은 SILVA138 분류기로 다시 붙인 것이 있으면 그걸 쓴다 (분류기 통일)
own_re = EXT / "own_taxonomy_silva138.tsv"
if own_re.exists():
    asv = pd.read_csv(ROOT / "data" / "own_asv_table.tsv", sep="\t", index_col=0) \
        if (ROOT / "data" / "own_asv_table.tsv").exists() else None
    if asv is not None:
        own_counts = collapse_to_genus(asv, pd.read_csv(own_re, sep="\t", index_col=0))
        W("  우리 데이터도 외부와 같은 SILVA 138 분류기로 재분류해 사용")
own_counts = own_counts[md.index]

# ---------------------------------------------------------------- 외부 시료 설계
meta = []
for s in ext_g.columns:
    m = re.match(r'^(SD\d+)\.(LS|NLS)\.(\d+)$', s)
    if m:
        meta.append((s, m.group(1), m.group(2), int(m.group(3))))
E = pd.DataFrame(meta, columns=['sample', 'subject', 'site', 'rep']).set_index('sample')
ext_g = ext_g[E.index]

W("=" * 78)
W("S2-1(2부)  외부 코호트 검증 — PRJNA1465550")
W(f"      시료 {ext_g.shape[1]}개 / 피험자 {E.subject.nunique()}명 / "
  f"병변 LS {int((E.site=='LS').sum())} · 비병변 NLS {int((E.site=='NLS').sum())}")
W(f"      속 {ext_g.shape[0]}개, 시료당 reads 중앙 {int(ext_g.sum(axis=0).median()):,}")
W("=" * 78)

# 기술 반복은 사람×부위 단위로 합친다 (짝 검정의 단위 = 사람)
key = E.subject + "|" + E.site
pool = ext_g.T.groupby(key.values).sum().T
P = pd.DataFrame([k.split("|") for k in pool.columns], columns=['subject', 'site'], index=pool.columns)
paired = [s for s in P.subject.unique()
          if {'LS', 'NLS'} <= set(P.site[P.subject == s])]
W(f"\n  기술 반복 합산 후: {pool.shape[1]}개 단위, 병변·비병변이 모두 있는 피험자 {len(paired)}명 -> 짝 검정 대상")

# 깊이 맞춤 (rarefy)
depth = int(pool.sum(axis=0).min())
W(f"  공통 깊이 {depth:,} reads 로 rarefy")
rng = np.random.default_rng(SEED)


def rarefy(col, d):
    v = col.values.astype(int); tot = v.sum()
    if tot <= d:
        return v.astype(float)
    idx = np.repeat(np.arange(len(v)), v)
    pick = rng.choice(idx, size=d, replace=False)
    return np.bincount(pick, minlength=len(v)).astype(float)


pool_r = pd.DataFrame({c: rarefy(pool[c], depth) for c in pool.columns}, index=pool.index)

# ---------------------------------------------------------------- 지표 계산
def logratio(tab, a, b):
    return np.log((tab.loc[a] + 0.5) / (tab.loc[b] + 0.5))


def shannon(col):
    p = col[col > 0] / col.sum()
    return float(-(p * np.log(p)).sum())


for need in ['Staphylococcus', 'Lawsonella', 'Cutibacterium']:
    if need not in pool_r.index:
        pool_r.loc[need] = 0.0
        W(f"  ! {need} 가 외부 코호트에 없어 0으로 채움 — 결과 해석 주의")

M = pd.DataFrame({
    'lipid_SL':  logratio(pool_r, 'Staphylococcus', 'Lawsonella') * -1,   # 부호: 높을수록 고피지형
    'lipid_CL':  logratio(pool_r, 'Cutibacterium', 'Lawsonella') * -1,
    'shannon':   [shannon(pool_r[c]) for c in pool_r.columns],
    'dominance': (pool_r.max(axis=0) / pool_r.sum(axis=0) * 100).values,
    'observed':  (pool_r > 0).sum(axis=0).values,
}, index=pool_r.columns).join(P)

# ---------------------------------------------------------------- 우리 모델 적용
common = own_counts.index.intersection(pool_r.index)
W(f"\n  두 코호트 공통 속 {len(common)}개 (우리 {own_counts.shape[0]} / 외부 {pool_r.shape[0]})")


def clr_on(tab, idx):
    x = np.log(tab.loc[idx].astype(float) + 0.5)
    return (x - x.mean(axis=0)).T


Xo, Xe = clr_on(own_counts, common).values, clr_on(pool_r, common).values
for target, col in [('Sebum', 'score_lipid'), ('Moisture', 'score_hydration')]:
    m = make_pipeline(StandardScaler(),
                      ElasticNetCV(l1_ratio=[.3, .5, .7, .9, 1.0], alphas=30,
                                   cv=5, random_state=SEED, max_iter=5000)).fit(Xo, md[target].values)
    M[col] = m.predict(Xe)

# ---------------------------------------------------------------- 짝 검정
W("\n### 병변(LS) vs 비병변(NLS) — 같은 사람 안에서 짝지은 검정 (Wilcoxon signed-rank)")
W(f"  {'지표':<26}{'NLS 중앙':>10}{'LS 중앙':>10}{'차이 중앙':>11}{'p':>9}   해석")
labels = {
    'lipid_SL':        ('지질 축  -log(Staph/Lawsonella)', '가설: LS에서 높다'),
    'lipid_CL':        ('지질 축  -log(Cuti/Lawsonella)', '가설: LS에서 높다'),
    'score_lipid':     ('지질 축 점수 (우리 모델)', '가설: LS에서 높다'),
    'shannon':         ('수분 축  Shannon 다양성', '음성대조: 차이 없어야'),
    'observed':        ('수분 축  관측 속 수', '음성대조: 차이 없어야'),
    'dominance':       ('수분 축  최우점 비율', '음성대조: 차이 없어야'),
    'score_hydration': ('수분 축 점수 (우리 모델)', '음성대조: 차이 없어야'),
}
res = {}
for k, (lab, hyp) in labels.items():
    a = M[(M.site == 'NLS') & M.subject.isin(paired)].set_index('subject')[k]
    b = M[(M.site == 'LS') & M.subject.isin(paired)].set_index('subject')[k]
    a, b = a.align(b, join='inner')
    d = b - a
    st = stats.wilcoxon(b, a)
    res[k] = (float(np.median(a)), float(np.median(b)), float(np.median(d)), float(st.pvalue))
    W(f"  {lab:<26}{np.median(a):>10.3f}{np.median(b):>10.3f}{np.median(d):>+11.3f}{st.pvalue:>9.4f}   {hyp}")

W("\n### 판정")
lip = [k for k in ['lipid_SL', 'lipid_CL', 'score_lipid'] if res[k][3] < 0.05 and res[k][2] > 0]
hyd = [k for k in ['shannon', 'observed', 'dominance', 'score_hydration'] if res[k][3] < 0.05]
W(f"  지질 축 지표 중 가설대로(LS↑, p<0.05): {len(lip)}/3  {lip}")
W(f"  수분 축 지표 중 유의하게 움직인 것: {len(hyd)}/4  {hyd if hyd else '없음 (음성대조 통과)'}")
if lip and not hyd:
    W("  -> 축 분리가 독립 코호트에서 재현됐다. 지질 축만 병변에서 움직이고 수분 축은 움직이지 않는다.")
elif lip and hyd:
    W("  -> 지질 축은 재현됐으나 수분 축도 함께 움직였다. 축 분리 주장은 약화된다 — 그대로 보고할 것.")
else:
    W("  -> 지질 축이 재현되지 않았다. 주장을 철회하고 그 사실을 포스터에 적을 것.")

M.to_csv(OUT / "05_external_metrics.tsv", sep="\t")
txt = "\n".join(rep)
(OUT / "05_external_validation.txt").write_text(txt, encoding="utf-8")
print(txt)
