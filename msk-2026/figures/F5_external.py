#!/usr/bin/env python
"""
F5  외부 검증 (2026-10-13 옛 F5 삭제로 F6 → F5) — 사전 예측은 실패했고, 그 실패가 두 축 구조를 설명한다.
 (A) 성적표 — 미리 건 예측(● 움직인다 / ○ 아니다)과 실제 결과. 여섯 항목 전부 어긋났다
 (B) 급내상관 — 지질 균형은 사람의 형질, 다양성은 부위의 상태
 (C) 지질 균형 log(Laws/Staph) ~ Shannon — 코호트별 ρ + 95% CI + 무선효과 통합

⚠️ v4(2026-10-11) 에서 (C) 를 바로잡았다.
   v3 의 'Ours ρ = 0.090, P = 0.481' 은 어느 산출물에도 없는 값이었다(P 는 다른 검정의 값과 같았다).
   v3 는 두 번째 코호트를 반대 부호 방향(log Staph/Laws)으로 그려 첫 코호트와 부호가 섞여 있었다.
   다시 계산하면 세 코호트 모두 양(+)이고 통합 ρ ≈ 0.20 [−0.04, 0.41] — "직교 재현"이 아니다.
   수치는 손으로 옮기지 않고 output/17_readout.tsv 에서 읽는다.
해석 문장은 그림에 없다 — 포스터 캡션에 있다.
수치 출처: output/07_reinterpret.txt · 10_baseline_only.txt · 17_readout.tsv
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (COL_W, CM, MOIST, SEBUM, INK, MUTED, GRID, RULE, FS_STAT, FS_ANNOT, FS_TICK,
                       ital, table, frame, head, fmt_p, num, ci_point, save_fig, save_panels)

RO = table("17_readout.tsv")
LS = f"log({ital('Lawsonella')}/{ital('Staphylococcus')})"
LC = f"log({ital('Lawsonella')}/{ital('Cutibacterium')})"


def scorecard(a):
    items = [LS, LC, "sebum model score", "Shannon", "observed genera", "hydration model score"]
    lipid = [True, True, True, False, False, False]
    pred = [True, True, True, False, False, False]
    # v8: baseline 짝 30명의 짝 검정 (output/19_baseline_paired.tsv). v7 의 p 는 06 의 baseline+treatment 합산(34명) 값이었다.
    T19 = table("19_baseline_paired.tsv")
    T19 = T19[T19.stratum == "baseline (.1)"].set_index('metric')
    pv = [float(T19.loc[k, 'p']) for k in
          ('lipid_SL', 'lipid_CL', 'score_lipid', 'shannon', 'observed', 'score_hydration')]
    obs = [p < 0.05 for p in pv]   # 채움 = 짝 검정 p < 0.05
    y = np.arange(len(items))[::-1]
    for yi, lp, pr, ob, p in zip(y, lipid, pred, obs, pv):
        cl = SEBUM if lp else MOIST
        a.plot([0, 1], [yi, yi], color=GRID, lw=3.0, zorder=2, solid_capstyle='round')
        for x0, flag in ((0, pr), (1, ob)):
            a.plot([x0], [yi], marker='o', ms=15, mfc=cl if flag else 'white', mec=cl, mew=2.4,
                   zorder=3)
        a.text(1.45, yi, fmt_p(p), ha='left', va='center', fontsize=FS_STAT,
               color=INK if p < .05 else MUTED)
    a.set_xlim(-0.40, 2.75); a.set_ylim(-0.6, len(items) - 0.4)
    a.set_xticks([0, 1, 1.98]); a.set_xticklabels(["predicted", "observed", "paired test"])
    a.set_yticks(y); a.set_yticklabels(items)
    for t, lp in zip(a.get_yticklabels(), lipid):
        t.set_color(SEBUM if lp else MOIST)
    a.tick_params(length=0)
    for s in a.spines.values():
        s.set_visible(False)


def icc(b):
    frame(b, grid='x')
    lab = [LS, LC, "Shannon", "observed genera"]
    val = [0.707, 0.788, 0.322, 0.444]
    cl = [SEBUM, SEBUM, MOIST, MOIST]
    y = np.arange(len(lab))[::-1]
    b.axvline(0.5, color=MUTED, lw=1.4, ls=(0, (4, 3)), zorder=1)
    for yi, v, c in zip(y, val, cl):
        b.plot([0, v], [yi, yi], color=c, lw=3.0, zorder=3, solid_capstyle='butt')
        b.plot([v], [yi], marker='o', ms=15, mfc=c, mec=c, zorder=4)
        b.text(v + 0.05, yi, f"{v:.2f}", ha='left', va='center', fontsize=FS_STAT, color=c,
               fontweight='bold')
    b.set_yticks(y); b.set_yticklabels(lab)
    for t, c in zip(b.get_yticklabels(), cl):
        t.set_color(c)
    b.set_xlim(0, 1.05); b.set_ylim(-0.6, len(lab) - 0.4)
    b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    b.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
    b.set_xlabel("ICC (paired sites)", labelpad=8)
    # v8b: site-/person-specific 해설은 포스터 범례로
    b.tick_params(axis='y', length=0)


def readout(c):
    frame(c, grid='x')
    spec = [("Ours", "unadjusted", "Ours (n = 27), unadjusted"),
            ("Ours", "age-adjusted", "Ours (n = 27), age-adjusted"),
            ("Cohort 1", "baseline pairs", "Cohort 1 (n = 30), baseline pairs"),
            ("Cohort 2", "published full-depth Shannon", "Cohort 2 (n = 21), full-depth Shannon"),
            ("Cohort 2", "Shannon of the other 28 genera", "Cohort 2 (n = 65), other 27 genera"),   # 표 머리 30속+Others 에서 4개 제외 = 27 (tsv spec 이름은 옛 표기 그대로)
            ("Pooled", "random effects", "Pooled, random effects")]
    y = np.array([6.0, 5.0, 4.0, 3.0, 2.0, 0.6])
    c.axvline(0, color=INK, lw=1.3, zorder=2)
    labels = []
    for yi, (co, sp, lab) in zip(y, spec):
        r = RO[(RO.cohort == co) & (RO.spec == sp)].iloc[0]
        in_pool = r.pooled >= 1
        col = INK if in_pool else MUTED
        if co == 'Pooled':
            c.fill([r.lo, r.rho, r.hi, r.rho], [yi, yi + 0.32, yi, yi - 0.32], color=INK, zorder=5)
        else:
            ci_point(c, yi, r.rho, r.lo, r.hi, col, r.p < 0.05, marker='s', ms=13)
        c.text(1.03, yi, f"{num(r.rho)}  [{num(r.lo, 2)}, {num(r.hi, 2)}]",
               transform=c.get_yaxis_transform(), ha='left', va='center', fontsize=FS_STAT,
               color=col, fontweight='bold' if co == 'Pooled' else 'normal')
        labels.append(lab)
    c.axhline(1.3, color=RULE, lw=1.0, zorder=1)
    c.set_yticks(y); c.set_yticklabels(labels)
    for t, (co, sp, _) in zip(c.get_yticklabels(), spec):
        r = RO[(RO.cohort == co) & (RO.spec == sp)].iloc[0]
        t.set_color(INK if r.pooled >= 1 else MUTED)
        if co == 'Pooled':
            t.set_fontweight('bold')
    c.set_xlim(-0.55, 0.92); c.set_ylim(0.0, 6.6)
    c.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
    c.set_xticklabels([num(v, 1) if v else "0" for v in (-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)])
    c.set_xlabel("Spearman $\\it{\\rho}$, 95% CI")
    c.tick_params(axis='y', length=0)


def build(name="F5_external_validation"):
    fig = plt.figure(figsize=(COL_W, 25 * CM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.62, left=0.30, right=0.985,
                             bottom=0.085, top=0.89)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 0.95], wspace=0.95)
    a = fig.add_subplot(top[0]); b = fig.add_subplot(top[1])
    c = fig.add_subplot(outer[1])
    pc = c.get_position()
    c.set_position([pc.x0, pc.y0, pc.width * 0.80, pc.height])   # 값 열 자리를 틀 밖에 남긴다
    scorecard(a); icc(b); readout(c)
    save_fig(fig, name)
    save_panels(fig, {'a': a, 'b': b, 'c': c}, "F5")
    return fig


if __name__ == "__main__":
    build()
