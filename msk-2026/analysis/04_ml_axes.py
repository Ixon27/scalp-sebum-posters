#!/usr/bin/env python
"""
S2-1(1부)  두 축 점수 — 군집 조성으로 피지·수분을 예측할 수 있는가.

정직성 장치 네 가지
 1. NS/SS 분류기는 만들지 않는다. 이 코호트에서 피지가 군을 100% 가르므로 순환논증이다.
    대신 피지를 '연속값'으로 예측한다.
 2. 그래도 군간 차이만 맞히는 것일 수 있으므로, 군 안에서만 본 상관을 따로 낸다.
    NS 안(피지 9–47)·SS 안(157–244)에서도 예측이 따라가면 진짜 용량-반응이다.
 3. 성능은 반드시 순열 귀무분포와 함께 보고한다. n=27에서 교차검증 R2 단독은 못 믿는다.
 4. 해석 가능한 기준선 — 로그비 log(Staphylococcus/Lawsonella) 하나만 쓴 모델과 비교한다.
    복잡한 모델이 이걸 못 이기면 로그비를 쓰는 게 맞다.
"""
import warnings
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

SEED = 20260825
NPERM = 499
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

g = pd.read_csv(OUT / "01_genus_counts.tsv", sep="\t", index_col=0)
md = pd.read_csv(OUT / "01_metadata.tsv", sep="\t", index_col=0)
g = g[md.index]

# 특징: 3개 이상 시료에서 검출된 속의 CLR (희소 잡음 제거, 그래도 넉넉히 남긴다)
feat = g.index[(g > 0).sum(axis=1) >= 3]
lg = np.log(g.loc[feat].astype(float) + 0.5)
X_all = (lg - lg.mean(axis=0)).T.values                  # 시료 x 속
names = list(feat)

lr = np.log((g.loc['Staphylococcus'] + 0.5) / (g.loc['Lawsonella'] + 0.5)).values.reshape(-1, 1)
grp = md.Group.values
rep = []; W = rep.append


def loo_predict(X, y, kind):
    """하나-빼기 교차검증 예측. 표준화·초매개변수 선택은 전부 학습 fold 안에서만."""
    pred = np.empty(len(y))
    for tr, te in LeaveOneOut().split(X):
        if kind == 'enet':
            m = make_pipeline(StandardScaler(),
                              ElasticNetCV(l1_ratio=[.3, .5, .7, .9, 1.0], alphas=30,
                                           cv=5, random_state=SEED, max_iter=5000))
        elif kind == 'rf':
            m = RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                      random_state=SEED, n_jobs=1)
        else:
            m = LinearRegression()
        m.fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return pred


def score(y, pred):
    ss_res = np.sum((y - pred) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return dict(R2=1 - ss_res / ss_tot, rho=stats.spearmanr(y, pred).statistic)


def perm_null(X, y, kind, n=NPERM):
    rng = np.random.default_rng(SEED)
    ys = [rng.permutation(y) for _ in range(n)]
    out = Parallel(n_jobs=-1)(delayed(lambda yy: score(yy, loo_predict(X, yy, kind))['R2'])(yy) for yy in ys)
    return np.array(out)


W("=" * 78)
W("S2-1(1부)  두 축 점수 — 조성으로 피지·수분을 예측할 수 있는가")
W(f"      특징 = 3시료 이상 검출된 {len(names)}속의 CLR / 하나-빼기 교차검증 / 순열 {NPERM}회")
W("=" * 78)

for target in ['Sebum', 'Moisture']:
    y = md[target].values.astype(float)
    W(f"\n{'='*70}\n### 목표: {target}\n{'='*70}")
    preds = {}
    for kind, lab, X in [('enet', 'Elastic Net (전체 조성)', X_all),
                         ('rf',   'Random Forest (전체 조성)', X_all),
                         ('lin',  'log(Staph/Lawsonella) 하나만', lr)]:
        p = loo_predict(X, y, kind); preds[lab] = p
        s = score(y, p)
        null = perm_null(X, y, kind)
        pv = (np.sum(null >= s['R2']) + 1) / (len(null) + 1)
        W(f"\n  {lab}")
        W(f"    교차검증 R2 = {s['R2']:+.3f}   Spearman(예측,실측) = {s['rho']:+.3f}")
        W(f"    순열 귀무 R2 중앙 {np.median(null):+.3f} (95% 상한 {np.percentile(null,95):+.3f})  ->  p = {pv:.4f}")
        # 군 안에서만 본 상관 — 군간 차이만 맞히는 게 아닌지 확인
        inner = []
        for gg in ['NS', 'SS']:
            m = grp == gg
            if m.sum() >= 5:
                r, pp = stats.spearmanr(y[m], p[m])
                inner.append(f"{gg}(n={m.sum()}) rho={r:+.2f} p={pp:.3f}")
        W(f"    군 안에서만: {' | '.join(inner)}")
    pd.DataFrame(preds, index=md.index).assign(observed=y, Group=grp)\
      .to_csv(OUT / f"04_pred_{target.lower()}.tsv", sep="\t")

# ---- 두 축 점수는 서로 독립인가
W(f"\n{'='*70}\n### 두 축 점수의 독립성\n{'='*70}")
ps = pd.read_csv(OUT / "04_pred_sebum.tsv", sep="\t", index_col=0)
pm = pd.read_csv(OUT / "04_pred_moisture.tsv", sep="\t", index_col=0)
col = 'Elastic Net (전체 조성)'
r, p = stats.spearmanr(ps[col], pm[col])
W(f"  피지 축 점수 ~ 수분 축 점수  rho = {r:+.3f}  p = {p:.4f}")

# ---- 어느 속이 각 축을 끌고 가나 (전체 데이터 적합 계수)
W(f"\n{'='*70}\n### 축별 기여 속 (전체 자료 적합 Elastic Net 계수, 상위 8개)\n{'='*70}")
ra_mean = (g / g.sum(axis=0) * 100).mean(axis=1)
for target in ['Sebum', 'Moisture']:
    y = md[target].values.astype(float)
    m = make_pipeline(StandardScaler(),
                      ElasticNetCV(l1_ratio=[.3, .5, .7, .9, 1.0], alphas=30,
                                   cv=5, random_state=SEED, max_iter=5000)).fit(X_all, y)
    co = pd.Series(m[-1].coef_, index=names)
    co = co[co != 0].reindex(co.abs().sort_values(ascending=False).index).dropna()
    W(f"\n  --- {target} (0이 아닌 계수 {len(co)}개 / {len(names)}속) ---")
    for k, v in co.head(8).items():
        W(f"    {k:<28} 계수 {v:+.3f}   (평균 상대존재량 {ra_mean.get(k,0):5.2f}%)")

txt = "\n".join(rep)
(OUT / "04_ml_axes.txt").write_text(txt, encoding="utf-8")
print(txt)
