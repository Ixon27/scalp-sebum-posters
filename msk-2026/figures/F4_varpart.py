#!/usr/bin/env python
"""
F4  분산분할 — 따로 재면 수분이 커 보이지만, 서로를 통제하면 피지가 가장 크다.
 (A) 기울기 그림 — 단독 R² 에서 고유 R² 로 갈 때 각 변수가 어디로 움직이나
 (B) 변수별 막대를 고유와 (나이와의) 공유로 쪼갠다 (adjusted R²)
v4(2026-10-11): KSABC 규격으로 다시 그렸다. 8-27 에 검토한 다른 형식(오일러·전체분산 띠)은 v4 에서 코드를 뺐다
— 모양은 시안 그림으로만 남아 있다(v3 스크립트는 v4 로 덮어썼다).
해석 문장은 그림에 넣지 않는다 — 포스터 캡션으로 뺀다.
수치 출처: output/02_varpart.txt (Bray–Curtis, partial dbRDA, 순열검정)
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (COL_W, CM, MOIST, SEBUM, AGE, INK, MUTED, GRID, DARK, PALE, FS_STAT, FS_ANNOT, FS_LAB,
                       frame, head, fmt_p, save_fig, save_panels)

NAMES = ['Sebum', 'Hydration', 'Age']
COLS = [SEBUM, MOIST, AGE]
ALONE = [0.1323, 0.2948, 0.2317]
UNIQ = [0.1489, 0.1235, 0.0518]
P_UNIQ = [0.0002, 0.0061, 0.0517]
P_ALONE = [0.0152, 0.0005, 0.0007]   # 02_varpart.txt [단독] — v8: 단독 점도 유의성대로 채운다 (v7 은 전부 속 빔)
U_S, U_H, U_A, SH_HA, RESID = 0.1376, 0.1101, 0.0324, 0.1907, 0.5692


def slope(a):
    frame(a, grid='y')
    for nm, cl, v0, v1, pu, pa in zip(NAMES, COLS, ALONE, UNIQ, P_UNIQ, P_ALONE):
        a.plot([0, 1], [v0, v1], color=cl, lw=3.0, zorder=3, solid_capstyle='round')
        a.plot([0], [v0], marker='o', ms=14, mfc=cl if pa < .05 else 'white', mec=cl, mew=2.4, zorder=4)
        a.plot([1], [v1], marker='o', ms=15, mfc=cl if pu < .05 else 'white', mec=cl, mew=2.4,
               zorder=5)
        a.text(-0.11, v0, f"{v0:.3f}", ha='right', va='center', fontsize=FS_STAT, color=cl,
               fontweight='bold')
        a.text(1.11, v1, f"{nm}  {v1:.3f}", ha='left', va='center', fontsize=FS_STAT, color=cl,
               fontweight='bold')
        a.text(2.55, v1, fmt_p(pu), ha='left', va='center', fontsize=FS_STAT,
               color=cl if pu < .05 else MUTED)
    a.set_xlim(-0.72, 3.40); a.set_ylim(0, 0.33)
    a.set_xticks([0, 1])
    a.set_xticklabels(["measured\nalone", "unique\n(others adjusted)"])
    a.set_ylabel("$\\it{R}$² (Bray–Curtis)")
    a.set_yticks(np.arange(0, 0.31, 0.10))


def split_bar(b):
    frame(b, grid='x')
    rows = [("Sebum", U_S, 0.0, SEBUM), ("Hydration", U_H, SH_HA, MOIST), ("Age", U_A, SH_HA, AGE)]
    y = np.arange(3)[::-1]
    for yi, (nm, u, sd, col) in zip(y, rows):
        b.barh([yi], [u], height=0.55, color=col, edgecolor='none', zorder=3)
        if sd:
            b.barh([yi], [sd], left=u, height=0.55, color=col, alpha=0.20, edgecolor=col, lw=1.2,
                   hatch='///', zorder=3)
            b.text(u + sd / 2, yi, f"{sd:.3f}", ha='center', va='center',
                   fontsize=FS_STAT, color=INK, zorder=5)
        if u >= 0.06:
            b.text(u / 2, yi, f"{u:.3f}", ha='center', va='center', fontsize=FS_STAT,
                   color='white', fontweight='bold', zorder=5)
        else:
            b.text(u / 2, yi + 0.30, f"{u:.3f}", ha='center', va='bottom', fontsize=FS_STAT,
                   color=col, fontweight='bold', zorder=5)
    b.set_yticks(y); b.set_yticklabels([r[0] for r in rows])
    for t, r in zip(b.get_yticklabels(), rows):
        t.set_color(r[3]); t.set_fontweight('bold'); t.set_fontsize(FS_LAB)
    b.set_xlim(0, 0.325); b.set_ylim(-0.6, 2.75)
    b.set_xlabel("adjusted $\\it{R}$² (Bray–Curtis)")
    b.tick_params(axis='y', length=0)


def build(name="F4_variance_partitioning"):
    fig = plt.figure(figsize=(COL_W, 14 * CM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.40, left=0.10, right=0.985,
                          bottom=0.20, top=0.78)
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])
    slope(a); split_bar(b)
    save_fig(fig, name)
    save_panels(fig, {'a': a, 'b': b}, "F4")
    return fig


if __name__ == "__main__":
    build()
