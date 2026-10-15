#!/usr/bin/env python
"""
F6  근거 지도 (2026-10-13 옛 F5 삭제로 F7 → F6) — 핵심 주장 5개를 세 가지 추정으로 다시 잰다 (2026-10-11 신설).
 행 = 주장, 줄 = 추정
   ● 전체 27명, 나이 통제 (정본)
   ■ 군 안 — 나이 + 군 통제, 순열은 군 안에서만. "피지 효과가 진단(NS/SS) 차이와 구별되는가"
   ◆ 사람 단위 15 — 같은 나이·같은 군을 한 사람으로 묶었을 때. "같은 사람이 여러 번 들어갔다면"
 선 = 95% CI (Fisher z + Bonett–Wright), 채움 = p < 0.05.
왜 이 그림인가: n = 27 에서 p 값 하나로는 "없다"와 "못 잡았다"를 가를 수 없고, 몇 점이 만든 결과인지도 안 보인다.
             CI 와 강건성 추정을 한 장에 두면 각 주장이 어디까지 버티는지가 바로 읽힌다.
해석 문장은 그림에 없다 — 포스터 캡션에 있다.
수치 출처: output/17_evidence.tsv (analysis/17_small_n/17_small_n.py)
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (COL_W, CM, MOIST, SEBUM, INK, MUTED, BAND, DARK, FS_STAT, FS_ANNOT, FS_LAB, FS_TICK,
                       RATIO_S, table, frame, head, num, ci_point, save_fig, save_panels)

EV = table("17_evidence_v8.tsv")   # v8: person 행 = 20_person_asv.py (ASV Shannon 시료값의 묶음 평균)
CLAIMS = [('indep', "Sebum vs hydration", INK),
          ('h_div', "Hydration vs Shannon", MOIST),
          ('s_div', "Sebum vs Shannon", SEBUM),
          ('s_bal', f"Sebum vs\n{RATIO_S}", SEBUM),
          ('h_bal', f"Hydration vs\n{RATIO_S}", MOIST)]
EST = [('all', 'o', "all"), ('within', 's', "within groups"), ('person', 'D', "age × group")]


def build(name="F6_evidence_map"):
    fig = plt.figure(figsize=(COL_W, 21 * CM))
    ax = fig.add_axes([0.375, 0.13, 0.43, 0.70])
    frame(ax, grid='x')
    gap = 1.0
    ys, ylabs = [], []
    for k, (key, lab, col) in enumerate(CLAIMS):
        base = -(k * (3 + gap))
        if k % 2 == 1:
            ax.axhspan(base - 2.55, base + 0.55, color=BAND, zorder=0, lw=0)
        for j, (est, mk, elab) in enumerate(EST):
            yi = base - j
            r = EV[(EV.claim == key) & (EV.estimate == est)].iloc[0]
            ci_point(ax, yi, r.rho, r.lo, r.hi, col, r.p < 0.05, marker=mk,
                     ms=13 if mk != 'D' else 11.5)
            ax.text(1.03, yi, f"{num(r.rho)}  [{num(r.lo, 2)}, {num(r.hi, 2)}]",
                    transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=FS_STAT,
                    color=col if r.p < 0.05 else MUTED,
                    fontweight='bold' if r.p < 0.05 else 'normal')
            ys.append(yi); ylabs.append(f"{elab} (n = {int(r.n)})")
        ax.text(-0.36, base - 1, lab, transform=ax.get_yaxis_transform(), ha='right', va='center',
                fontsize=FS_LAB, fontweight='bold', color=col)
    ax.axvline(0, color=INK, lw=1.3, zorder=2)
    ax.set_yticks(ys); ax.set_yticklabels(ylabs, fontsize=FS_TICK, color=INK)
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(min(ys) - 0.7, 0.7)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xticks([-0.8, -0.4, 0, 0.4, 0.8])
    ax.set_xticklabels([num(v, 1) if v else "0" for v in (-0.8, -0.4, 0, 0.4, 0.8)])
    ax.set_xlabel("partial $\\it{\\rho}$ (age-adjusted), 95% CI")
    ax.annotate("partial $\\it{\\rho}$  [95% CI]", xy=(1.03, 1), xycoords='axes fraction',
                xytext=(0, 6), textcoords='offset points', ha='left', va='bottom',
                fontsize=FS_STAT, color=MUTED, annotation_clip=False)
    hs = [Line2D([], [], marker=mk, ls='-', lw=2.4, ms=12 if mk != 'D' else 10.5, color=DARK,
                 mfc=DARK, label=t)
          for _, mk, t in (('all', 'o', "all 27 participants"),
                           ('within', 's', "within groups (NS/SS adjusted)"),
                           ('person', 'D', "age × group clusters (15 units)"))]
    # v8b: 기호 key 는 포스터 범례 문장으로 (hs 는 남겨 둠)
    save_fig(fig, name)
    save_panels(fig, {'': ax}, "F6")
    return fig


if __name__ == "__main__":
    build()
