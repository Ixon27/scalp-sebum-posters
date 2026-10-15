#!/usr/bin/env python
"""
서버 전용 통합 분석 — (1) 두 축 예측 모델 + 순열검정  (2) 독립 코호트 외부 검증

로컬 Mac에서 돌리지 말 것. 순열 999회 x 하나-빼기 27겹 x 내부 5겹 CV 는 맥을 멈춘다.
실행:
  conda activate qiime2-amplicon-2024.5
  cd <서버의 msk-2026 폴더>
  setsid nohup python analysis/06_server_ml_external.py > output/06.log 2>&1 < /dev/null & disown

정직성 장치
 1. NS/SS 분류기는 만들지 않는다 — 이 코호트에서 피지가 군을 100% 가르므로 순환논증.
 2. 군 안에서만 본 상관을 따로 낸다 — 군간 차이만 맞히는 게 아닌지 확인.
 3. 성능은 순열 귀무분포와 함께 보고한다.
 4. 해석 가능한 기준선 log(Staphylococcus/Lawsonella) 하나와 비교한다.
 5. 외부 검증은 음성대조(수분 축)를 함께 건다 — 둘 다 움직이면 축 분리 주장이 틀린 것.
"""
import os, re, warnings, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed

warnings.filterwarnings("ignore")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

SEED, NPERM, NJOBS = 20260825, 999, 32
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"; OUT.mkdir(exist_ok=True)
rep = []


def W(*a):
    s = "".join(str(x) for x in a)
    rep.append(s); print(s, flush=True)


import inspect
_ALPHA_KW = 'n_alphas' if 'n_alphas' in inspect.signature(ElasticNetCV).parameters else 'alphas'


def enet():
    kw = {_ALPHA_KW: 20}
    return make_pipeline(StandardScaler(),
                         ElasticNetCV(l1_ratio=[.5, .9, 1.0], cv=5,
                                      max_iter=2000, tol=1e-3, selection='random',
                                      random_state=SEED, n_jobs=1, **kw))


def make_model(kind):
    if kind == 'enet':
        return enet()
    if kind == 'rf':
        return RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                     random_state=SEED, n_jobs=1)
    return LinearRegression()


def loo_predict(X, y, kind):
    pred = np.empty(len(y))
    for tr, te in LeaveOneOut().split(X):
        m = make_model(kind); m.fit(X[tr], y[tr]); pred[te] = m.predict(X[te])
    return pred


def r2(y, p):
    return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)


def rank_genus(taxon):
    for p in (x.strip() for x in taxon.split(';')):
        if p.startswith('g__'):
            return p[3:]
    return ''


def to_genus(counts, tax):
    gg = tax['Taxon'].map(rank_genus)
    gg = gg[gg != '']
    idx = counts.index.intersection(gg.index)
    return counts.loc[idx].join(gg.rename('genus')).groupby('genus').sum()


# =================================================================== 자료
tax_own = pd.read_csv(ROOT / "data/external/own_taxonomy_silva138.tsv", sep="\t", index_col=0)
asv_own = pd.read_csv(ROOT / "data/own_asv_table.tsv", sep="\t", index_col=0)
md = pd.read_csv(ROOT / "output/01_metadata.tsv", sep="\t", index_col=0)
own = to_genus(asv_own, tax_own)[md.index]

ext_counts = pd.read_csv(ROOT / "data/external/feature-table.tsv", sep="\t", skiprows=1, index_col=0)
ext_tax = pd.read_csv(ROOT / "data/external/taxonomy.tsv", sep="\t", index_col=0)
ext = to_genus(ext_counts, ext_tax)

W("=" * 78)
W("서버 통합 분석 — 두 축 예측 모델 + 독립 코호트 외부 검증")
W(f"  우리 코호트 : {own.shape[1]}시료 x {own.shape[0]}속 (NS {int((md.Group=='NS').sum())} / SS {int((md.Group=='SS').sum())})")
W(f"  외부 코호트 : {ext.shape[1]}시료 x {ext.shape[0]}속  (PRJNA1465550)")
W(f"  순열 {NPERM}회 / 병렬 {NJOBS} / 난수고정 {SEED}")
W("=" * 78)

feat = own.index[(own > 0).sum(axis=1) >= 3]
lg = np.log(own.loc[feat].astype(float) + 0.5)
X_all = (lg - lg.mean(axis=0)).T.values
lr_own = np.log((own.loc['Staphylococcus'] + 0.5) / (own.loc['Lawsonella'] + 0.5)).values.reshape(-1, 1)
grp = md.Group.values

# =================================================================== (1) 예측 모델
W("\n" + "#" * 70)
W("# (1) 조성으로 피지·수분을 예측할 수 있는가")
W("#" * 70)
W(f"  특징 = 3시료 이상 검출된 {len(feat)}속의 CLR")

for target in ['Sebum', 'Moisture']:
    y = md[target].values.astype(float)
    W(f"\n{'='*66}\n### 목표: {target}\n{'='*66}")
    for kind, lab, X in [('enet', 'Elastic Net (전체 조성)', X_all),
                         ('rf', 'Random Forest (전체 조성)', X_all),
                         ('lin', 'log(Staph/Lawsonella) 하나만', lr_own)]:
        p = loo_predict(X, y, kind)
        obs = r2(y, p); rho = stats.spearmanr(y, p).statistic
        rng = np.random.default_rng(SEED)
        ys = [rng.permutation(y) for _ in range(NPERM)]
        null = np.array(Parallel(n_jobs=NJOBS)(
            delayed(lambda yy: r2(yy, loo_predict(X, yy, kind)))(yy) for yy in ys))
        pv = (np.sum(null >= obs) + 1) / (NPERM + 1)
        W(f"\n  {lab}")
        W(f"    교차검증 R2 = {obs:+.3f}   Spearman(예측,실측) = {rho:+.3f}")
        W(f"    순열 귀무 R2 중앙 {np.median(null):+.3f}, 95% 상한 {np.percentile(null,95):+.3f}  ->  p = {pv:.4f}")
        inner = []
        for gg in ['NS', 'SS']:
            m = grp == gg
            r, pp = stats.spearmanr(y[m], p[m])
            inner.append(f"{gg}(n={m.sum()}) rho={r:+.2f} p={pp:.3f}")
        W(f"    군 안에서만: {' | '.join(inner)}   <- 군간 차이만 맞히는 게 아닌지 확인")
        pd.DataFrame({'obs': y, 'pred': p, 'Group': grp}, index=md.index)\
          .to_csv(OUT / f"06_pred_{target.lower()}_{kind}.tsv", sep="\t")

# 기여 속
W(f"\n{'='*66}\n### 축별 기여 속 (전체 자료 적합 Elastic Net, 0이 아닌 계수)\n{'='*66}")
ra_own = (own / own.sum(axis=0) * 100).mean(axis=1)
for target in ['Sebum', 'Moisture']:
    m = enet().fit(X_all, md[target].values.astype(float))
    co = pd.Series(m[-1].coef_, index=list(feat))
    co = co[co != 0]
    co = co.reindex(co.abs().sort_values(ascending=False).index)
    W(f"\n  --- {target} ({len(co)}속 선택 / {len(feat)}속 중) ---")
    for k, v in co.head(10).items():
        W(f"    {k:<30} 계수 {v:+8.3f}   평균 상대존재량 {ra_own.get(k,0):6.2f}%")

# =================================================================== (2) 외부 검증
W("\n" + "#" * 70)
W("# (2) 독립 코호트 외부 검증 — PRJNA1465550 (지루성 피부염 환자 병변/비병변 짝)")
W("#" * 70)

meta = [(s, m.group(1), m.group(2)) for s in ext.columns
        if (m := re.match(r'^(SD\d+)\.(LS|NLS)\.\d+$', s))]
E = pd.DataFrame(meta, columns=['sample', 'subject', 'site']).set_index('sample')
ext = ext[E.index]
W(f"  시료 {ext.shape[1]} / 피험자 {E.subject.nunique()}명 / 병변 {int((E.site=='LS').sum())} · 비병변 {int((E.site=='NLS').sum())}")

pool = ext.T.groupby((E.subject + "|" + E.site).values).sum().T
P = pd.DataFrame([c.split("|") for c in pool.columns], columns=['subject', 'site'], index=pool.columns)
paired = sorted(s for s in P.subject.unique() if {'LS', 'NLS'} <= set(P.site[P.subject == s]))
W(f"  기술반복 합산 후 {pool.shape[1]}단위, 병변·비병변 모두 있는 피험자 {len(paired)}명")

depth = int(pool.sum(axis=0).min())
rng = np.random.default_rng(SEED)


def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d:
        return v.astype(float)
    pick = rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False)
    return np.bincount(pick, minlength=len(v)).astype(float)


pool_r = pd.DataFrame({c: rarefy(pool[c].values, depth) for c in pool.columns}, index=pool.index)
W(f"  공통 깊이 {depth:,} reads 로 rarefy")

for need in ['Staphylococcus', 'Cutibacterium', 'Lawsonella']:
    if need not in pool_r.index:
        pool_r.loc[need] = 0.0
        W(f"  ! 경고: {need} 가 외부 코호트에 없다")
    else:
        nz = int((pool_r.loc[need] > 0).sum())
        W(f"  {need:<16} 검출 {nz}/{pool_r.shape[1]}단위, 평균 상대존재량 "
          f"{(pool_r.loc[need]/pool_r.sum(axis=0)*100).mean():5.2f}%")


def shannon(v):
    p = v[v > 0] / v.sum()
    return float(-(p * np.log(p)).sum())


M = pd.DataFrame({
    'lipid_SL': -np.log((pool_r.loc['Staphylococcus'] + .5) / (pool_r.loc['Lawsonella'] + .5)),
    'lipid_CL': -np.log((pool_r.loc['Cutibacterium'] + .5) / (pool_r.loc['Lawsonella'] + .5)),
    'shannon': [shannon(pool_r[c].values) for c in pool_r.columns],
    'observed': (pool_r > 0).sum(axis=0).astype(float),
    'dominance': pool_r.max(axis=0) / pool_r.sum(axis=0) * 100,
}, index=pool_r.columns).join(P)

common = sorted(set(own.index) & set(pool_r.index))
W(f"\n  두 코호트 공통 속 {len(common)}개 (우리 {own.shape[0]} / 외부 {pool_r.shape[0]})")


def clr_of(tab):
    x = np.log(tab.loc[common].astype(float) + 0.5)
    return (x - x.mean(axis=0)).T.values


Xo, Xe = clr_of(own), clr_of(pool_r)
for target, col in [('Sebum', 'score_lipid'), ('Moisture', 'score_hydration')]:
    m = enet().fit(Xo, md[target].values.astype(float))
    M[col] = m.predict(Xe)
    W(f"  {target} 모델(공통 속 재학습) 자기 코호트 교차검증 R2 = "
      f"{r2(md[target].values.astype(float), loo_predict(Xo, md[target].values.astype(float), 'enet')):+.3f}")

W("\n### 병변(LS) vs 비병변(NLS) — 같은 사람 안에서 짝지은 Wilcoxon signed-rank")
W(f"  {'지표':<34}{'NLS':>9}{'LS':>9}{'차이':>9}{'p':>9}   가설")
LAB = {
    'lipid_SL': ('지질축 -log(Staph/Lawsonella)', 'LS에서 높다'),
    'lipid_CL': ('지질축 -log(Cuti/Lawsonella)', 'LS에서 높다'),
    'score_lipid': ('지질축 점수 (우리 모델)', 'LS에서 높다'),
    'shannon': ('수분축 Shannon 다양성', '음성대조: 차이 없어야'),
    'observed': ('수분축 관측 속 수', '음성대조: 차이 없어야'),
    'dominance': ('수분축 최우점 비율', '음성대조: 차이 없어야'),
    'score_hydration': ('수분축 점수 (우리 모델)', '음성대조: 차이 없어야'),
}
res = {}
for k, (lab, hyp) in LAB.items():
    a = M[(M.site == 'NLS') & M.subject.isin(paired)].set_index('subject')[k].sort_index()
    b = M[(M.site == 'LS') & M.subject.isin(paired)].set_index('subject')[k].sort_index()
    a, b = a.align(b, join='inner')
    st = stats.wilcoxon(b.values, a.values)
    res[k] = (np.median(a), np.median(b), np.median(b.values - a.values), st.pvalue)
    W(f"  {lab:<34}{np.median(a):>9.3f}{np.median(b):>9.3f}"
      f"{np.median(b.values-a.values):>+9.3f}{st.pvalue:>9.4f}   {hyp}")

lip = [k for k in ['lipid_SL', 'lipid_CL', 'score_lipid'] if res[k][3] < 0.05 and res[k][2] > 0]
hyd = [k for k in ['shannon', 'observed', 'dominance', 'score_hydration'] if res[k][3] < 0.05]
W("\n### 판정")
W(f"  지질 축 지표 중 가설대로 (LS 높음, p<0.05): {len(lip)}/3 {lip}")
W(f"  수분 축 지표 중 유의하게 움직인 것: {len(hyd)}/4 {hyd if hyd else '없음 — 음성대조 통과'}")
if lip and not hyd:
    W("  -> 축 분리가 독립 코호트에서 재현됐다.")
elif lip and hyd:
    W("  -> 지질 축은 재현됐으나 수분 축도 움직였다. 축 분리 주장 약화 — 그대로 보고.")
else:
    W("  -> 지질 축이 재현되지 않았다. 주장을 철회하고 그대로 적을 것.")

M.to_csv(OUT / "06_external_metrics.tsv", sep="\t")
(OUT / "06_server_ml_external.txt").write_text("\n".join(rep), encoding="utf-8")
print("\n저장 완료: output/06_server_ml_external.txt")
