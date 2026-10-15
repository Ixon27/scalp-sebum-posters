#!/usr/bin/env python
"""
F1  두 축 평면 — 포스터의 뼈대.
가로 = 피지 실측(µg/cm²), 세로 = 수분(나이 통제 순위 잔차). 27명을 평면에 찍고 양쪽에 군별 밀도를 붙인다.
가로로는 두 군이 완전히 갈리고(빈 구간 47→157), 세로로는 겹친다.
해석 문장은 그림에 없다 — 포스터 캡션에 있다.
수치 출처: output/01_two_axis.txt [1] (피지 ~ 수분 나이 통제 −0.015, p = 0.945)
"""
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
os.environ['FS_SCALE'] = '1.4'   # v8d: 포스터 배율 약 0.72 → 실제 글자 약 15 pt
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (COL_W, CM, MOIST, SEBUM, INK, MUTED, NS_COL, SS_COL, NS_T, SS_T, FS_STAT,
                       FS_ANNOT, load, resid, frame, group_scatter, group_legend, stat, fmt_p,
                       save_fig, save_panels)

md, alpha, lr, counts = load()
sebum = md['Sebum'].to_numpy(float)
mr = resid(md['Moisture'], md['Age'])
grp = md['Group'].to_numpy()
ns, ss = grp == 'NS', grp == 'SS'
p_grp = stats.mannwhitneyu(mr[ns], mr[ss]).pvalue   # 수분 잔차의 군 차이 (이 패널에서만 쓴다)


def build(name="F1_two_axes"):
    fig = plt.figure(figsize=(COL_W, 15 * CM))
    gs = fig.add_gridspec(2, 2, width_ratios=[5, 1], height_ratios=[0.8, 5],
                          wspace=0.03, hspace=0.03, left=0.10, right=0.98, bottom=0.13, top=0.91)
    axm = fig.add_subplot(gs[1, 0])
    axt = fig.add_subplot(gs[0, 0], sharex=axm)
    axr = fig.add_subplot(gs[1, 1], sharey=axm)

    XLIM, YLIM = (-14, 264), (-19.5, 19.5)
    lo, hi = sebum[ss].min(), sebum[ss].max()

    # ━━ 본 패널
    frame(axm)
    axm.axvspan(lo - 8, hi + 8, color=SS_T, alpha=0.35, zorder=0, lw=0)
    axm.axhline(0, color=MUTED, lw=1.2, ls=(0, (5, 4)), zorder=1)
    group_scatter(axm, sebum, mr, grp, s=230)
    axm.set_xlim(*XLIM); axm.set_ylim(*YLIM)
    axm.set_yticks([-15, -10, -5, 0, 5, 10, 15])
    axm.set_xlabel("Sebum (µg/cm²)", color=SEBUM, labelpad=8)
    axm.set_ylabel("Hydration\n(age-adjusted rank residual)", color=MOIST, labelpad=6)   # v8d: 두 줄 — 축 높이보다 길다
    # 통계 — 데이터가 비어 있는 왼쪽 위 (KSABC 처럼 패널 안)
    axm.text(0.025, 0.965, f"Sebum vs hydration: {stat(rho=-0.015, p=0.945)}\n"
                           f"Hydration, NS vs SS: {fmt_p(p_grp)}",
             transform=axm.transAxes, ha='left', va='top', fontsize=FS_STAT, color=INK,
             linespacing=1.6, zorder=8)

    # ━━ 위 = 피지 밀도 (두 군이 완전히 갈린다)
    gx = np.linspace(*XLIM, 500); top = 0
    for sel, c, t in ((ns, NS_COL, NS_T), (ss, SS_COL, SS_T)):
        d = stats.gaussian_kde(sebum[sel], bw_method=0.45)(gx)
        axt.fill_between(gx, d, color=t, alpha=0.6, lw=0, zorder=2)
        axt.plot(gx, np.where(d > d.max() * 0.01, d, np.nan), color=c, lw=2.2, zorder=3)
        top = max(top, d.max())
    axt.set_ylim(0, top * 1.15)
    # v8b: 'no overlap' 괄호·문구는 포스터 범례로 옮김

    # ━━ 오른쪽 = 수분 밀도 (두 군이 겹친다)
    gy = np.linspace(*YLIM, 500); right = 0
    for sel, c, t in ((ns, NS_COL, NS_T), (ss, SS_COL, SS_T)):
        d = stats.gaussian_kde(mr[sel], bw_method=0.55)(gy)
        axr.fill_betweenx(gy, d, color=t, alpha=0.6, lw=0, zorder=2)
        axr.plot(np.where(d > d.max() * 0.01, d, np.nan), gy, color=c, lw=2.2, zorder=3)
        right = max(right, d.max())
    axr.set_xlim(0, right * 1.2)

    for ax in (axt, axr):
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    # 위 밀도의 rug — 축 바로 위 짧은 눈금
    for sel, c in ((ns, NS_COL), (ss, SS_COL)):
        axt.vlines(sebum[sel], 0, top * 0.16, color=c, lw=1.8, zorder=5)
        axr.hlines(mr[sel], 0, right * 0.12, color=c, lw=1.8, zorder=5)

    group_legend(fig, x=0.98, y=1.0)
    save_fig(fig, name)
    save_panels(fig, {'': (axm, axt, axr)}, "F1")
    return fig


if __name__ == "__main__":
    build()
