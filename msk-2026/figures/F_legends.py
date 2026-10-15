#!/usr/bin/env python
"""
포스터용 범례 조각 (v8, 2026-10-14) — 패널 이미지에는 figure 범례가 빠지므로 따로 저장해 포스터에 붙인다.
 groups : NS / SS (세로 2줄)
 F3     : NS / SS / Lawsonella not detected + Sebum / Hydration (세로)
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (NS_COL, SS_COL, EDGE, DARK, FS_ANNOT, FIGDIR, DATE, VER, ital, axis_handles)

G = [Line2D([], [], marker='o', ls='', ms=13, mfc=NS_COL, mec=EDGE, mew=0.9, label="NS"),
     Line2D([], [], marker='o', ls='', ms=13, mfc=SS_COL, mec=EDGE, mew=0.9, label="SS")]
ND = Line2D([], [], marker='o', ls='', ms=13, mfc='white', mec=NS_COL, mew=2.0,
            label=f"{ital('Lawsonella')}\nnot detected")


def save(handles, tag, ncol=1, title=None):
    fig = plt.figure(figsize=(3, 3))
    lg = fig.legend(handles=handles, loc='center', ncol=ncol, columnspacing=1.6, frameon=False, fontsize=FS_ANNOT,
                    handletextpad=0.4, labelspacing=0.7, borderpad=0, title=title,
                    title_fontsize=FS_ANNOT, alignment='left')
    fig.canvas.draw()
    bb = lg.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).padded(0.05)
    for ext in ('pdf', 'png'):
        fig.savefig(f"{FIGDIR}/panels/{DATE}_{tag}_v{VER}.{ext}", bbox_inches=bb,
                    dpi=300 if ext == 'png' else None)
    plt.close(fig)


save(G, "legend_groups")
save(axis_handles(sebum_first=True), "legend_F3")   # v8d: NS·SS 는 F3A 띠 라벨 색이 key — D 용 Sebum·Hydration 만
# v8d: F6 모양 key — 추정 수준 (원·사각형·마름모). 색 key 와 같은 취급으로 그림 옆에 둔다
F6K = [Line2D([], [], marker=mk, ls='-', lw=2.4, ms=12 if mk != 'D' else 10.5, color=DARK,
              mfc=DARK, label=t)
       for mk, t in (('o', "all 27 participants"), ('s', "within groups (NS/SS adjusted)"),
                     ('D', "age × group clusters (15 units)"))]
save(F6K, "legend_F6", ncol=3)
# v8e: F5·F6 모양 key — 그림 오른쪽 빈 칸에 세로로 (두 줄 이름)
from fig_style import MUTED, INK
def mk(m, filled=True, col=DARK, ls='', label='', ms=13):
    return Line2D([], [], marker=m, ls=ls, lw=2.4, ms=ms, color=col, mfc=col if filled else 'white',
                  mec=col, mew=2.2, label=label)
save([mk('o', True, label="change expected\nor $\\it{p}$ < 0.05"),
      mk('o', False, label="no change expected\nor not significant")], "legend_F5a",
     title="(A)")   # v8f: B 옆에 놓이므로 A 용임을 밝힌다
save([mk('s', True, ls='-', label="$\\it{p}$ < 0.05"),
      mk('s', False, ls='-', label="$\\it{p}$ ≥ 0.05"),
      mk('s', False, col=MUTED, ls='-', label="sensitivity,\nnot pooled"),
      mk('D', True, col=INK, label="pooled estimate", ms=15)], "legend_F5c")
save([mk('o', ls='-', label="all 27\nparticipants"),
      mk('s', ls='-', label="within groups\n(NS/SS adjusted)"),
      mk('D', ls='-', label="age × group\nclusters (15 units)", ms=11.5)], "legend_F6v")
print("저장: panels/legend_groups · legend_F3 · legend_F6 · legend_F5a · legend_F5c · legend_F6v")
