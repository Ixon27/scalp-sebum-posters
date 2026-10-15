#!/usr/bin/env python
"""
F2  수분 축 — 수분이 다양성을 지배하고, 피지도 군도 그 자리를 건드리지 못한다.
 (A) Shannon ~ 수분 (나이 통제 순위 잔차) + 부트스트랩 95% 신뢰대
 (B) 다양성 5지표 × 두 축 — 편상관 + 95% CI (2026-10-11 v4: 덤벨 → CI 로 교체)
 (C) 군 비교 — 지루성군의 다양성이 정상과 다르지 않다 (KSABC F1a 처럼 박스 + 점 + 반쪽 밀도)
 (D) 수분과 홍반을 서로 통제했을 때 무엇이 살아남나 (정직성 항목)
해석 문장은 그림에 넣지 않는다 — 포스터 캡션으로 뺐다.
수치 출처: output/01_two_axis.txt [2] · 17_ci_alpha.tsv · 16_canonical_redness.txt B
"""
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import os
os.environ['FS_SCALE'] = '1.5'   # v8d: 포스터 배율 0.63–0.70 → 실제 글자 약 15 pt
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (COL_W, CM, MOIST, MOIST_T, SEBUM, INK, MUTED, OTHER, FRAME, EDGE, NS_COL, SS_COL, NS_T,
                       SS_T, FS_STAT, FS_ANNOT, FS_TICK, load, resid, table, frame, head, inside_top,
                       group_scatter, group_legend, fit_band, ci_point, axis_handles, stat, fmt_p,
                       num, save_fig, save_panels)

md, alpha, lr, counts = load()
grp = md['Group'].to_numpy(); ns, ss = grp == 'NS', grp == 'SS'
mr = resid(md['Moisture'], md['Age'])
sh = alpha['Shannon'].to_numpy(float)
CI = table("17_ci_alpha.tsv")
IDX = ['Shannon', 'Simpson', 'ACE', 'Observed', 'Dominance']
LAB = ['Shannon', 'Simpson', 'ACE', 'Observed\nASVs', 'Dominance']   # 표시용 — 지표는 ASV 수준(01_two_axis.py)


def build(name="F2_hydration_axis"):
    fig = plt.figure(figsize=(COL_W, 30 * CM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.08], height_ratios=[1, 1], wspace=0.42,
                          hspace=0.62, left=0.09, right=0.985, bottom=0.10, top=0.89)
    a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    c, d = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    # ━━━━━ (A) 수분 ~ Shannon
    frame(a)
    fit_band(a, mr, sh, MOIST, MOIST_T)
    group_scatter(a, mr, sh, grp)
    a.set_xlabel("Hydration (age-adjusted rank residual)", color=MOIST)
    a.set_ylabel("Shannon index")
    a.set_ylim(1.80, 4.45)
    head(a, None, None, stat(rho=0.553, p=0.003))

    # ━━━━━ (B) 5지표 × 두 축, 95% CI
    frame(b, grid='x')
    y0 = np.arange(len(IDX))[::-1].astype(float)
    b.axvline(0, color=INK, lw=1.3, zorder=2)
    for yi, idx in zip(y0, IDX):
        for axis, col, dy in (('Hydration', MOIST, 0.23), ('Sebum', SEBUM, -0.23)):   # v8d: 글자가 커져 간격 넓힘
            r = CI[(CI['index'] == idx) & (CI['axis'] == axis)].iloc[0]
            ci_point(b, yi + dy, r.rho, r.lo, r.hi, col, r.p < 0.05, ms=12)
            b.text(0.985, yi + dy, num(r.rho), transform=b.get_yaxis_transform(),
                   ha='right', va='center', fontsize=FS_ANNOT,
                   color=col if r.p < 0.05 else MUTED,
                   fontweight='bold' if r.p < 0.05 else 'normal')
    b.set_yticks(y0); b.set_yticklabels(LAB)
    b.set_xlim(-1.0, 1.22); b.set_ylim(-0.6, len(IDX) - 0.4)
    b.set_xticks([-0.8, -0.4, 0, 0.4, 0.8])
    b.set_xticklabels([num(v, 1) if v else "0" for v in (-0.8, -0.4, 0, 0.4, 0.8)])
    b.tick_params(axis='y', length=0)
    b.set_xlabel("partial $\\it{\\rho}$ (age-adjusted), 95% CI")
    b.legend(handles=axis_handles(), loc='upper center', bbox_to_anchor=(0.42, -0.20), ncol=2,
             handletextpad=0.4, columnspacing=1.4, borderpad=0, fontsize=FS_ANNOT)

    # ━━━━━ (C) 군 비교 — KSABC F1a 형식
    frame(c, grid='y')
    rng = np.random.default_rng(20260825)
    for x0, sel, col, tint in ((0, ns, NS_COL, NS_T), (1, ss, SS_COL, SS_T)):
        v = sh[sel]
        q1, q2, q3 = np.percentile(v, [25, 50, 75])
        iqr = q3 - q1
        wlo, whi = v[v >= q1 - 1.5 * iqr].min(), v[v <= q3 + 1.5 * iqr].max()
        # 반쪽 밀도 (오른쪽)
        # 밀도는 관측 범위 안에서만 그린다 — 범위 밖으로 늘이면 없는 값(Shannon < 2)을 그리게 된다
        gy = np.linspace(v.min(), v.max(), 200)
        kd = stats.gaussian_kde(v, bw_method=0.6)(gy)
        c.fill_betweenx(gy, x0 + 0.08, x0 + 0.08 + kd / kd.max() * 0.34, color=tint, alpha=0.9,
                        lw=0, zorder=2)
        c.plot(x0 + 0.08 + kd / kd.max() * 0.34, gy, color=col, lw=2.2, zorder=2)
        # 박스
        c.plot([x0, x0], [wlo, whi], color=FRAME, lw=1.6, zorder=3)
        c.add_patch(plt.Rectangle((x0 - 0.06, q1), 0.12, iqr, facecolor=col, edgecolor=FRAME,
                                  lw=1.4, zorder=4))
        c.plot([x0 - 0.06, x0 + 0.06], [q2, q2], color='white', lw=2.0, zorder=5)
        c.scatter([x0], [v.mean()], marker='D', s=110, color=INK, zorder=6)
        # 점 (왼쪽)
        c.scatter(x0 - 0.24 + rng.normal(0, 0.035, sel.sum()), v, s=120, color=col,
                  edgecolor=EDGE, lw=0.9, zorder=3)
    c.set_xticks([0, 1]); c.set_xticklabels(['NS', 'SS'])
    c.set_xlim(-0.55, 1.6); c.set_ylim(1.80, 4.75)
    c.set_ylabel("Shannon index")
    inside_top(c, fmt_p(0.246), col=MUTED)

    # ━━━━━ (D) 수분 vs 홍반, 서로 통제  (출처: 16_canonical_redness.txt B)
    cols = ['Hydration\n(erythema-adjusted)', 'Erythema\n(hydration-adjusted)']
    val = np.array([[+0.547, +0.320], [+0.549, +0.359], [+0.573, +0.552],
                    [+0.508, +0.483], [-0.554, -0.356]])
    pv = np.array([[0.0025, 0.1051], [0.0039, 0.0672], [0.0018, 0.0022],
                   [0.0060, 0.0101], [0.0026, 0.0710]])
    base = [MOIST, OTHER]
    n = len(IDX)
    for i in range(n):
        for j in range(2):
            inten = min(abs(val[i, j]) / 0.60, 1.0) * 0.60
            sig = pv[i, j] < 0.05
            # 채움 = 유의, 속 빔(흰 면 + 회색 테두리) = 비유의 — CI 점과 같은 규칙 (v5)
            d.add_patch(plt.Rectangle((j + 0.06, n - 0.94 - i), 0.88, 0.88,
                                      facecolor=to_rgba(base[j], max(inten, 0.30)) if sig else 'white',
                                      edgecolor=base[j] if sig else MUTED, lw=2.0 if sig else 1.4,
                                      zorder=2))
            d.text(j + 0.5, n - 0.5 - i, num(val[i, j]), ha='center', va='center',
                   fontsize=FS_STAT, color=INK if sig else MUTED,
                   fontweight='bold' if sig else 'normal', zorder=3)
    d.set_xlim(0, 2); d.set_ylim(0, n)
    d.set_xticks([0.5, 1.5]); d.set_xticklabels(cols, fontsize=FS_TICK)
    d.set_yticks(np.arange(n) + 0.5); d.set_yticklabels(LAB[::-1])
    d.tick_params(length=0)
    for s in d.spines.values():
        s.set_visible(False)
    d.set_xlabel("mutually adjusted partial $\\it{\\rho}$",
                 labelpad=10)

    group_legend(fig, x=0.985, y=1.0)
    save_fig(fig, name)
    save_panels(fig, {'a': a, 'b': b, 'c': c, 'd': d}, "F2")
    return fig


if __name__ == "__main__":
    build()
