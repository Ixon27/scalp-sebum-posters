#!/usr/bin/env python
"""
17  작은 표본(n=27) 보강 — p 값 대신 효과크기 + 95% 신뢰구간, 그리고 강건성 검정.

왜: 핵심 주장 중 셋이 "없다"는 주장(두 축 독립 · 피지↛다양성 · 수분↛균형)이다.
    n=27 에서 p>0.05 는 "없다"가 아니라 "못 잡았다"이므로, 얼마나 큰 효과까지 배제되는지를 CI 로 보인다.
    또 피지는 정상(9–47)과 지루성(157–244)을 겹침 0 으로 가르므로, 피지 효과가 군 차이와 구별되는지 따로 본다.

계산 (전부 Spearman 순위 기반, 나이 통제는 01 과 같은 순위 잔차 방식)
  [A] 다양성 5지표 × 두 축 편상관 + 95% CI                    -> 17_ci_alpha.tsv   (F2B)
  [B] 검출 기준 4속 × 두 축 편상관 + 95% CI (ρ·q 는 01 산출물)  -> 17_ci_genus.tsv   (F3D)
  [C] 핵심 주장 5개 × 세 가지 추정 (전체 · 군 안 · 사람 단위)    -> 17_evidence.tsv   (F7)
  [D] 한 명씩 빼기 · 나이 묶음 하나씩 빼기 (수분 축 · 피지 축)     -> 보고서에만
  [E] 지질 균형 ~ 다양성, 코호트 3개 + 무선효과 통합              -> 17_readout.tsv    (F6C)

CI: Fisher z 변환 + Bonett & Wright(2000) Spearman 분산  se = sqrt((1 + ρ²/2) / (n − 3 − k)),  k = 통제 변수 수.
    공식 대조: Ruscio 2008 식 7 · Bishara & Hittner 2017 식 16 (원문 비공개). 편상관 확장의 포함률은 18_ci_coverage.py 로 확인.
p : 순열 9,999회, 난수 20260825. 나이 통제는 01 과 같이 x 를 섞는다.
    군 통제(군 안) 검정은 x 를 **군 안에서만** 섞는다 — 피지는 군과 완전히 얽혀 있어 전체를 섞으면 귀무가 틀린다.
부호: 지질 균형은 외부 코호트 스크립트(07·10)와 같은 log(Lawsonella/Staphylococcus) 방향으로 통일한다.
      08·09(두 번째 코호트)는 log(Staphylococcus/Lawsonella) 방향이라 여기서 부호를 뒤집어 맞춘다.

실행: python3 analysis/17_small_n/17_small_n.py   (로컬 수 분, 병렬 없음)
"""
import re
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
OUT = f"{ROOT}/output"
SEED, NPERM = 20260825, 9999
RNG = np.random.default_rng(SEED)
sys.path.insert(0, f"{ROOT}/analysis/08_ext2_PRJNA1199187")
from xlsx_loader import load as load_xlsx   # noqa: E402

LOG = []


def say(s=""):
    print(s)
    LOG.append(s)


# ━━━━━━━━━━ 통계 도구 ━━━━━━━━━━
def rk(v):
    return stats.rankdata(v)


def pcor(x, y, covs=()):
    """편 Spearman — 순위를 통제변수 순위(군은 0/1 그대로)에 회귀한 잔차끼리의 상관."""
    rx, ry = rk(x), rk(y)
    Z = np.column_stack([np.ones(len(rx))] + [c if set(np.unique(c)) <= {0.0, 1.0} else rk(c)
                                              for c in covs])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def perm_p(x, y, covs=(), strata=None):
    obs = pcor(x, y, covs)
    null = np.empty(NPERM)
    for i in range(NPERM):
        if strata is None:
            xp = RNG.permutation(x)
        else:
            xp = np.array(x, float).copy()
            for s in np.unique(strata):
                m = strata == s
                xp[m] = RNG.permutation(xp[m])
        null[i] = pcor(xp, y, covs)
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)


def ci(r, n, k=0):
    se = np.sqrt((1 + r * r / 2) / (n - 3 - k))
    z = np.arctanh(r)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def dl_pool(rs, ns, ks):
    """DerSimonian–Laird 무선효과 통합 (Fisher z, Bonett–Wright 분산)."""
    rs, ns, ks = map(np.asarray, (rs, ns, ks))
    z = np.arctanh(rs)
    v = (1 + rs ** 2 / 2) / (ns - 3 - ks)
    w = 1 / v
    zf = np.sum(w * z) / w.sum()
    Q = float(np.sum(w * (z - zf) ** 2)); df = len(rs) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - np.sum(w ** 2) / w.sum()))
    ws = 1 / (v + tau2)
    zr = np.sum(ws * z) / ws.sum(); se = 1 / np.sqrt(ws.sum())
    p = 2 * stats.norm.sf(abs(zr / se))
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return dict(rho=float(np.tanh(zr)), lo=float(np.tanh(zr - 1.96 * se)),
                hi=float(np.tanh(zr + 1.96 * se)), p=float(p), Q=Q, Qp=float(stats.chi2.sf(Q, df)),
                I2=I2, tau2=tau2)


def fmt(r, lo, hi, p=None):
    s = f"ρ = {r:+.3f}  [{lo:+.2f}, {hi:+.2f}]"
    return s + (f"  p = {p:.4f}" if p is not None else "")


# ━━━━━━━━━━ 자료 ━━━━━━━━━━
md = pd.read_csv(f"{OUT}/01_metadata.tsv", sep="\t", index_col=0)
al = pd.read_csv(f"{OUT}/01_alpha.tsv", sep="\t", index_col=0).loc[md.index]
lr = pd.read_csv(f"{OUT}/01_logratios.tsv", sep="\t", index_col=0).loc[md.index]
cnt = pd.read_csv(f"{OUT}/01_genus_counts.tsv", sep="\t", index_col=0)[md.index]
N = len(md)
AGE = md['Age'].to_numpy(float)
SEB = md['Sebum'].to_numpy(float)
MOI = md['Moisture'].to_numpy(float)
GRP = (md['Group'] == 'SS').to_numpy(float)
SH = al['Shannon'].to_numpy(float)
SL = lr['log(Staphylococcus/Lawsonella)'].to_numpy(float)   # 01 방향 (Staph 위)
LS = -SL                                                     # 지질 균형 방향 (Lawsonella 위)

say("=" * 78)
say("17  작은 표본 보강 — 효과크기 + 95% CI · 군 안 추정 · 사람 단위 · 영향점 · 코호트 통합")
say(f"    n = {N} (NS {int((GRP == 0).sum())} / SS {int(GRP.sum())}) · 순열 {NPERM}회 · 난수 {SEED}")
say("    CI = Fisher z + Bonett–Wright Spearman 분산, se = sqrt((1+ρ²/2)/(n−3−k))")
say("=" * 78)

# ━━━━━━━━━━ [A] 다양성 5지표 × 두 축 ━━━━━━━━━━
say("\n[A] 다양성 5지표 × 두 축 — 나이 통제 편상관 + 95% CI  (-> 17_ci_alpha.tsv, F2B)")
rowsA = []
for idx in ['Shannon', 'Simpson', 'ACE', 'Observed', 'Dominance']:
    y = al[idx].to_numpy(float)
    for axis, x in (('Hydration', MOI), ('Sebum', SEB)):
        r, p = perm_p(x, y, (AGE,))
        lo, hi = ci(r, N, 1)
        rowsA.append(dict(index=idx, axis=axis, rho=r, lo=lo, hi=hi, p=p, n=N))
        say(f"    {idx:10s} {axis:9s} {fmt(r, lo, hi, p)}")
pd.DataFrame(rowsA).to_csv(f"{OUT}/17_ci_alpha.tsv", sep="\t", index=False, float_format="%.6f")

# 01 정본과 대조 — 반올림 2자리 값이 어긋나면 멈춘다
ref01 = {('Shannon', 'Hydration'): 0.55, ('Simpson', 'Hydration'): 0.55, ('ACE', 'Hydration'): 0.55,
         ('Observed', 'Hydration'): 0.50, ('Dominance', 'Hydration'): -0.56,
         ('Shannon', 'Sebum'): 0.29, ('Simpson', 'Sebum'): 0.36, ('ACE', 'Sebum'): 0.13,
         ('Observed', 'Sebum'): 0.18, ('Dominance', 'Sebum'): -0.37}
for r in rowsA:
    assert abs(round(r['rho'], 2) - ref01[(r['index'], r['axis'])]) <= 0.011, r

# ━━━━━━━━━━ [B] 속 단위 ━━━━━━━━━━
say("\n[B] 검출 기준 4속 × 두 축 — ρ·q 는 01 산출물, CI 만 여기서  (-> 17_ci_genus.tsv, F3D)")
rowsB = []
for axis, f in (('Sebum', '01_genus_sebum.tsv'), ('Hydration', '01_genus_moisture.tsv')):
    g = pd.read_csv(f"{OUT}/{f}", sep="\t", index_col=0)
    for gen, row in g.iterrows():
        lo, hi = ci(row['rho'], N, 1)
        rowsB.append(dict(genus=gen, axis=axis, rho=row['rho'], lo=lo, hi=hi, p=row['p'],
                          q=row['q'], n=N))
        say(f"    {gen:15s} {axis:9s} {fmt(row['rho'], lo, hi, row['p'])}  q = {row['q']:.4f}")
pd.DataFrame(rowsB).to_csv(f"{OUT}/17_ci_genus.tsv", sep="\t", index=False, float_format="%.6f")

# ━━━━━━━━━━ [C] 핵심 주장 5개 × 세 가지 추정 ━━━━━━━━━━
say("\n[C] 핵심 주장 5개 × 세 가지 추정  (-> 17_evidence.tsv, F7)")
say("    전체     = 27명, 나이 통제 (정본)")
say("    군 안    = 27명, 나이 + 군 통제, 순열은 군 안에서만  ← 피지 효과가 군 차이와 구별되는가")
say("    사람 단위 = 같은 나이·같은 군을 한 사람으로 묶은 15단위, 나이 통제 (13_conservative.txt (b))")

# 사람 단위 값은 13 산출물에서 읽는다 — 손으로 옮겨 적지 않는다
txt13 = open(f"{OUT}/13_conservative.txt", encoding="utf-8").read()
sec_b = txt13.split("(b) 사람 단위")[1].split("(c) 사람 단위")[0]


def from13(label):
    m = re.search(re.escape(label) + r"\s*\(나이 통제\)\s+rho=([+-][\d.]+)\s+p=([\d.]+)", sec_b)
    if m is None:
        raise SystemExit(f"13_conservative.txt (b) 에서 '{label}' 를 못 찾았다")
    return float(m.group(1)), float(m.group(2))


CLAIMS = [  # (key, 라벨, x, y, 13 라벨, 축)
    ('indep', 'Sebum vs hydration', SEB, MOI, '피지 ~ 수분', 'both'),
    ('h_div', 'Hydration → Shannon', MOI, SH, '수분 ~ Shannon', 'Hydration'),
    ('s_div', 'Sebum → Shannon', SEB, SH, '피지 ~ Shannon', 'Sebum'),
    ('s_bal', 'Sebum → log(Staph/Laws)', SEB, SL, '피지 ~ log(Staph/Laws)', 'Sebum'),
    ('h_bal', 'Hydration → log(Staph/Laws)', MOI, SL, '수분 ~ log(Staph/Laws)', 'Hydration'),
]
rowsC = []
for key, lab, x, y, l13, axis in CLAIMS:
    r1, p1 = perm_p(x, y, (AGE,))
    r2, p2 = perm_p(x, y, (AGE, GRP), strata=GRP)
    r3, p3 = from13(l13)
    for est, r, p, n, k in (('all', r1, p1, N, 1), ('within', r2, p2, N, 2),
                            ('person', r3, p3, 15, 1)):
        lo, hi = ci(r, n, k)
        rowsC.append(dict(claim=key, label=lab, axis=axis, estimate=est, rho=r, lo=lo, hi=hi,
                          p=p, n=n, k=k))
    say(f"\n    {lab}")
    for est in ('all', 'within', 'person'):
        rr = [z for z in rowsC if z['claim'] == key and z['estimate'] == est][0]
        say(f"      {est:7s} n={rr['n']:2d}  {fmt(rr['rho'], rr['lo'], rr['hi'], rr['p'])}")
pd.DataFrame(rowsC).to_csv(f"{OUT}/17_evidence.tsv", sep="\t", index=False, float_format="%.6f")

say("\n    군별로 따로 (참고 — 표본이 작아 CI 가 넓다)")
for g, m in (('NS', GRP == 0), ('SS', GRP == 1)):
    for lab, x, y in (('피지 ~ log(Staph/Laws)', SEB, SL), ('수분 ~ Shannon', MOI, SH),
                      ('피지 ~ Shannon', SEB, SH)):
        r = pcor(x[m], y[m], (AGE[m],)); lo, hi = ci(r, int(m.sum()), 1)
        say(f"      {g} (n={int(m.sum())})  {lab:24s} {fmt(r, lo, hi)}")

law0 = (cnt.loc['Lawsonella'] == 0).to_numpy()
say(f"\n    Lawsonella 미검출(0 reads): NS {int(law0[GRP == 0].sum())}/16 · SS {int(law0[GRP == 1].sum())}/11")
say("    -> 미검출 시료의 log(Staph/Laws) 는 pseudocount(0.5)가 정한 값이다. 정상군 대부분이 여기에 해당한다.")

# ━━━━━━━━━━ [D] 영향점 ━━━━━━━━━━
say("\n[D] 영향점 — 한 명씩 빼기 · 나이 묶음(같은 나이 2~4명) 하나씩 빼기")
for lab, x, y in (('수분 ~ Shannon (나이 통제)', MOI, SH), ('피지 ~ log(Staph/Laws) (나이 통제)', SEB, SL)):
    full = pcor(x, y, (AGE,))
    loo = np.array([pcor(np.delete(x, i), np.delete(y, i), (np.delete(AGE, i),)) for i in range(N)])
    say(f"\n    {lab}: 전체 {full:+.3f} · 한 명 빼기 범위 [{loo.min():+.3f}, {loo.max():+.3f}]"
        f" (최저: {md.index[int(loo.argmin())]} 제외)")
    worst = (None, 1.0, 0.0)
    for ag in sorted(np.unique(AGE)):
        keep = AGE != ag
        if (~keep).sum() < 2:
            continue
        r, p = perm_p(x[keep], y[keep], (AGE[keep],))
        say(f"      {int(ag)}세 묶음 제외 ({','.join(md.index[~keep])})  ρ = {r:+.3f}  p = {p:.4f}")
        if abs(r) < abs(worst[2]) or worst[0] is None:
            worst = (int(ag), p, r)
    say(f"      -> 가장 약해질 때: {worst[0]}세 묶음 제외, ρ = {worst[2]:+.3f}, p = {worst[1]:.4f}")

# ━━━━━━━━━━ [E] 지질 균형 ~ 다양성, 코호트 3개 ━━━━━━━━━━
say("\n[E] 지질 균형 log(Lawsonella/Staphylococcus) ~ Shannon — 코호트 3개  (-> 17_readout.tsv, F6C)")
say("    ⚠️ 2026-08-27 판 F6C 의 'Ours ρ = 0.090, P = 0.481' 은 어느 산출물에도 없는 값이었다.")
say("       P = 0.481 은 다른 검정(수분 ~ log(Staph/Laws), p = 0.4811)의 값과 같다.")
say("    ⚠️ 같은 판 F6C 는 두 번째 코호트 값을 log(Staph/Laws) 방향(−0.136)으로 그려 첫 코호트(log(Laws/Staph))와 부호가 섞여 있었다.")

rowsE = []
r, p = perm_p(LS, SH)
lo, hi = ci(r, N, 0)
rowsE.append(dict(cohort='Ours', spec='unadjusted', rho=r, lo=lo, hi=hi, p=p, n=N, k=0, pooled=1))
r, p = perm_p(LS, SH, (AGE,))
lo, hi = ci(r, N, 1)
rowsE.append(dict(cohort='Ours', spec='age-adjusted', rho=r, lo=lo, hi=hi, p=p, n=N, k=1, pooled=0))

# 첫 외부 코호트 — 10_baseline_only.txt 의 baseline 짝 30명 (lipid_SL 은 이미 log(Laws/Staph))
txt10 = open(f"{OUT}/10_baseline_only.txt", encoding="utf-8").read()
sec = txt10.split("### baseline만 (.1)")[1].split("### treatment만")[0]
m = re.search(r"lipid_SL ~ shannon\s+rho = ([+-][\d.]+)\s+p = ([\d.]+).*?\(n=(\d+)\)", sec)
r, p, n = float(m.group(1)), float(m.group(2)), int(m.group(3))
lo, hi = ci(r, n, 0)
rowsE.append(dict(cohort='Cohort 1', spec='baseline pairs', rho=r, lo=lo, hi=hi, p=p, n=n, k=0,
                  pooled=1))

# 두 번째 외부 코호트 — 부록표에서 다시 계산, 도시 통제, 부호를 log(Laws/Staph) 로 뒤집는다
d2 = load_xlsx(f"{ROOT}/data/PRJNA1199187/fmicb-2025-1595030-TableS1.xlsx")
g = d2['top30.sample.g']; hdr = g[0][1:]; ids = [row[0] for row in g[1:]]
B = np.array([[float(v) if v not in ('', None) else 0.0 for v in row[1:1 + len(hdr)]] for row in g[1:]])
meta = {row[0]: row for row in d2['metadata'][1:]}
city = np.array([1.0 if meta[s][3] == 'Shantou' else 0.0 for s in ids])
col = lambda nm: B[:, hdr.index(nm)]
pc = min(col('Staphylococcus')[col('Staphylococcus') > 0].min(),
         col('Lawsonella')[col('Lawsonella') > 0].min()) / 2
LS2 = np.log((col('Lawsonella') + pc) / (col('Staphylococcus') + pc))


def shannon(M):
    P = M / M.sum(1, keepdims=True)
    return np.array([-(q[q > 0] * np.log(q[q > 0])).sum() for q in P])


ab = d2['alpha_index bacteria']; apub = {row[0]: float(row[6]) for row in ab[1:]}
com = [s for s in ids if s in apub]; ix = [ids.index(s) for s in com]
r, p = perm_p(LS2[ix], np.array([apub[s] for s in com]), (city[ix],))
lo, hi = ci(r, len(ix), 1)
rowsE.append(dict(cohort='Cohort 2', spec='published full-depth Shannon', rho=r, lo=lo, hi=hi, p=p,
                  n=len(ix), k=1, pooled=1))
keep = [i for i, c in enumerate(hdr) if c not in ('Cutibacterium', 'Staphylococcus', 'Lawsonella',
                                                   'Others')]
R = B[:, keep]; ok = R.sum(1) > 0
r, p = perm_p(LS2[ok], shannon(R[ok]), (city[ok],))
lo, hi = ci(r, int(ok.sum()), 1)
rowsE.append(dict(cohort='Cohort 2', spec='Shannon of the other 28 genera', rho=r, lo=lo, hi=hi,
                  p=p, n=int(ok.sum()), k=1, pooled=0))

for z in rowsE:
    say(f"    {z['cohort']:9s} {z['spec']:32s} n={z['n']:2d}  {fmt(z['rho'], z['lo'], z['hi'], z['p'])}")

pl = [z for z in rowsE if z['pooled']]
P = dl_pool([z['rho'] for z in pl], [z['n'] for z in pl], [z['k'] for z in pl])
say(f"\n    통합 (무선효과, 같은 양끼리: 우리 비보정 + 코호트1 + 코호트2 전체 깊이)")
say(f"      {fmt(P['rho'], P['lo'], P['hi'], P['p'])}   Q = {P['Q']:.2f} (p = {P['Qp']:.3f}), I² = {P['I2']:.0%}")
rowsE.append(dict(cohort='Pooled', spec='random effects', rho=P['rho'], lo=P['lo'], hi=P['hi'],
                  p=P['p'], n=sum(z['n'] for z in pl), k=0, pooled=2))
alt = [rowsE[0], rowsE[2], rowsE[4]]
P2 = dl_pool([z['rho'] for z in alt], [z['n'] for z in alt], [z['k'] for z in alt])
say(f"    민감도 (코호트2 를 잔여 28속 Shannon n=65 로 바꾸면)")
say(f"      {fmt(P2['rho'], P2['lo'], P2['hi'], P2['p'])}   I² = {P2['I2']:.0%}")
ext = [rowsE[2], rowsE[3]]
P3 = dl_pool([z['rho'] for z in ext], [z['n'] for z in ext], [z['k'] for z in ext])
say(f"    외부 두 코호트만")
say(f"      {fmt(P3['rho'], P3['lo'], P3['hi'], P3['p'])}   I² = {P3['I2']:.0%}")
pd.DataFrame(rowsE).to_csv(f"{OUT}/17_readout.tsv", sep="\t", index=False, float_format="%.6f")

say("\n    해석: 세 코호트 모두 같은 방향(양)이고 통합 CI 가 0 근처까지 걸친다.")
say("          '직교가 재현됐다'가 아니라 '미생물 지표 둘(지질 균형·다양성)은 약하게 얽혀 있고 정밀도가 낮다'가 맞는 서술이다.")
say("          두 축의 독립은 생리 지표(피지 ~ 수분, 우리 코호트에서만 측정)에 대한 주장이고, 외부 코호트는 그것을 직접 검정하지 못한다.")

with open(f"{OUT}/17_small_n.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
print(f"\n저장: output/17_small_n.txt · 17_ci_alpha.tsv · 17_ci_genus.tsv · 17_evidence.tsv · 17_readout.tsv")
