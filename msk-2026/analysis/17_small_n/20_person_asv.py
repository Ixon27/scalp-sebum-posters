#!/usr/bin/env python
"""
20  F6 '나이×군 묶음' 추정을 전체·군 안 추정과 같은 지표 정의로 다시 (2026-10-14).

왜: 17_evidence.tsv 의 person 행은 13_conservative.txt (b) 에서 읽는다. 13 (b) 는 속 count 를 묶음별로
    평균한 조성을 다시 rarefy 해 Shannon 을 새로 계산한다(속 수준). 반면 all·within 행의 Shannon 은
    01_alpha.tsv(ASV 수준) 다. 그래서 F6 의 세 추정을 비교하면 '묶음 효과'와 '지표 정의 차이'가 섞인다 (Codex 검토).
무엇을: 같은 나이 + 같은 군 = 한 단위(15개, 13 과 같은 정의). 각 단위의 값은 그 단위 시료들의
    **시료별 지표 평균** — Shannon 은 01_alpha.tsv(ASV), 로그비는 01_logratios.tsv, 임상값은 01_metadata.
    all·within 과 같은 지표를 같은 정의로 쓰고, 바뀌는 것은 '묶음' 하나뿐이다.
    통계: 나이 통제 편 Spearman, 순열 9,999회, CI = Fisher z + Bonett–Wright (n = 15, k = 1) — 17 과 같다.
산출: output/20_person_asv.txt · output/17_evidence_v8.tsv (17_evidence 의 person 행만 이 값으로 바꾼 판, F6 v8 이 읽는다)
"""
import numpy as np, pandas as pd
from scipy import stats

ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
OUT = f"{ROOT}/output"
SEED, NPERM = 20260825, 9999
RNG = np.random.default_rng(SEED)
LOG = []


def say(s=""):
    print(s); LOG.append(s)


def pcor(x, y, covs=()):   # 17 과 같은 정의
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    Z = np.column_stack([np.ones(len(rx))] + [stats.rankdata(c) for c in covs])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def perm_p(x, y, covs=()):
    obs = pcor(x, y, covs)
    null = np.array([pcor(RNG.permutation(x), y, covs) for _ in range(NPERM)])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)


def ci(r, n, k=0):
    se = np.sqrt((1 + r * r / 2) / (n - 3 - k)); z = np.arctanh(r)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


md = pd.read_csv(f"{OUT}/01_metadata.tsv", sep="\t", index_col=0)
al = pd.read_csv(f"{OUT}/01_alpha.tsv", sep="\t", index_col=0).loc[md.index]
lr = pd.read_csv(f"{OUT}/01_logratios.tsv", sep="\t", index_col=0).loc[md.index]
D = pd.DataFrame({'Sebum': md['Sebum'], 'Moisture': md['Moisture'], 'Age': md['Age'],
                  'Shannon': al['Shannon'], 'SL': lr['log(Staphylococcus/Lawsonella)'],
                  'unit': md['Age'].astype(str) + "|" + md['Group']})
U = D.groupby('unit')[['Sebum', 'Moisture', 'Age', 'Shannon', 'SL']].mean()
n = len(U)
say("=" * 78)
say(f"20  나이×군 묶음 {n}단위 — 시료별 지표(ASV Shannon · 로그비 · 임상값)의 단위 평균, 나이 통제")
say("=" * 78)
CL = [('indep', 'Sebum', 'Moisture'), ('h_div', 'Moisture', 'Shannon'), ('s_div', 'Sebum', 'Shannon'),
      ('s_bal', 'Sebum', 'SL'), ('h_bal', 'Moisture', 'SL')]
ev = pd.read_csv(f"{OUT}/17_evidence.tsv", sep="\t")
for key, xv, yv in CL:
    r, p = perm_p(U[xv].to_numpy(float), U[yv].to_numpy(float), (U['Age'].to_numpy(float),))
    lo, hi = ci(r, n, 1)
    m = (ev.claim == key) & (ev.estimate == 'person')
    old = ev.loc[m].iloc[0]
    say(f"  {key:6s} {xv:>8} ~ {yv:<8}  ρ = {r:+.3f} [{lo:+.2f}, {hi:+.2f}]  p = {p:.4f}"
        f"     (13 (b) 판: ρ = {old.rho:+.3f}, p = {old.p:.4f})")
    ev.loc[m, ['rho', 'lo', 'hi', 'p', 'n']] = [r, lo, hi, p, n]
ev.to_csv(f"{OUT}/17_evidence_v8.tsv", sep="\t", index=False, float_format="%.6f")
open(f"{OUT}/20_person_asv.txt", "w", encoding="utf-8").write("\n".join(LOG))
print("저장: output/20_person_asv.txt · output/17_evidence_v8.tsv")
