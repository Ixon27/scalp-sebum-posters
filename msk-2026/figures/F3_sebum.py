#!/usr/bin/env python
"""
F3  피지 축 — 피지가 우점 3속의 균형을 지배한다.
 (A) 조성 바코드 — 27명을 피지 순으로 세워 놓고 3속 조성을 쌓는다. 아래 띠가 피지 실측값
 (B) log(Staph/Lawsonella) ~ 피지 실측   (C) 같은 로그비 ~ 수분 (음성대조, y축 공유)
 (D) 속 단위 — 두 축의 편상관 + 95% CI (2026-10-11 v4: 덤벨 → CI 로 교체)

v4 에서 바꾼 것과 이유
 - B 의 가로축을 순위 잔차에서 **피지 실측(µg/cm²)** 으로 바꿨다. 순위 잔차는 두 군 사이 빈 구간(47→157)을
   눌러 연속 기울기처럼 보이게 만든다. 실측축에서는 상관의 대부분이 군 사이 격차라는 것이 그대로 보인다
   (군 안 추정은 F6). 통계는 그대로 나이 통제 편상관이다.
 - 적합선을 뺐다. 데이터가 없는 빈 구간을 가로지르는 선은 그 사이에 기울기가 있다고 말해 버린다.
 - Lawsonella 가 0 reads 인 시료(정상 16명 중 10명)는 속 빈 점. 그 로그비는 pseudocount(0.5)가 정한 값이다.
해석 문장은 그림에 없다 — 포스터 캡션에 있다.
수치 출처: output/01_two_axis.txt [3]–[4] · 17_ci_genus.tsv
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
import os
os.environ['FS_SCALE'] = '1.5'   # v8d: 포스터 배율 0.63–0.70 → 실제 글자 약 15 pt
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from fig_style import (COL_W, CM, MOIST, SEBUM, INK, MUTED, RULE, FRAME, LW_FRAME, SEBUM_RAMP, NS_COL,
                       SS_COL, GEN_COL, FS_ANNOT,
                       FS_TICK, load, resid, table, frame, head, group_scatter, group_legend,
                       ci_point, axis_handles, stat, num, ital, RATIO, save_fig, save_panels)

md, alpha, lr, counts = load()
grp = md['Group'].to_numpy()
mr = resid(md['Moisture'], md['Age'])
SL = lr['log(Staphylococcus/Lawsonella)'].to_numpy(float)
LAW0 = (counts.loc['Lawsonella'] == 0).to_numpy()
GCI = table("17_ci_genus.tsv")


def build(name="F3_sebum_axis"):
    fig = plt.figure(figsize=(COL_W, 31 * CM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.05], hspace=0.42, left=0.085,
                             right=0.985, bottom=0.08, top=0.90)
    top = outer[0].subgridspec(2, 1, height_ratios=[6, 0.9], hspace=0.06)
    bot = outer[1].subgridspec(1, 2, width_ratios=[2.0, 0.92], wspace=0.42)
    bc = bot[0].subgridspec(1, 2, wspace=0.17)
    a = fig.add_subplot(top[0]); strip = fig.add_subplot(top[1], sharex=a)
    b = fig.add_subplot(bc[0]); c = fig.add_subplot(bc[1], sharey=b); d = fig.add_subplot(bot[1])

    # ━━━━━ (A) 조성 바코드 — 피지 순
    frame(a, grid='y')
    ra = counts / counts.sum(0) * 100
    order = md['Sebum'].sort_values().index
    three = ['Staphylococcus', 'Cutibacterium', 'Lawsonella']
    comp = ra.loc[three, order].T
    comp['other genera'] = 100 - comp.sum(1)
    x = np.arange(len(order)); bottom = np.zeros(len(order))
    for g in three + ['other genera']:
        a.bar(x, comp[g], bottom=bottom, width=0.84, color=GEN_COL[g], edgecolor='white', lw=0.6,
              zorder=3, label=ital(g) if g != 'other genera' else 'other genera')
        bottom += comp[g].to_numpy()
    nb = (md.loc[order, 'Group'] == 'NS').sum() - 0.5
    a.axvline(nb, color=INK, lw=2.0, zorder=5)
    a.set_xlim(-0.7, len(order) - 0.3); a.set_ylim(0, 100)
    a.set_ylabel("Relative\nabundance (%)")   # v8d: 두 줄 — 아래 Sebum 띠 이름과 겹침
    a.tick_params(axis='x', bottom=False, labelbottom=False)
    # v8: 포스터에서 F3A 는 0.73 배로 들어가 범례가 작다 → 속 범례만 크게 (22 pt → 포스터 약 16 pt)
    a.legend(loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=4, handlelength=1.2,
             handleheight=1.0, handletextpad=0.45, columnspacing=1.6, borderpad=0.2,
             fontsize=25)   # v8e: 포스터 배율 0.585 → 약 14.6 pt

    cmap = LinearSegmentedColormap.from_list('sebum', SEBUM_RAMP)
    sv = md.loc[order, 'Sebum'].to_numpy(float)
    strip.imshow(sv[None, :], aspect='auto', cmap=cmap, norm=Normalize(0, 250),
                 extent=(-0.7, len(order) - 0.3, 0, 1))
    strip.set_yticks([0.5]); strip.set_yticklabels(["Sebum"], color=SEBUM, fontweight='bold')
    strip.tick_params(length=0, labelbottom=False)
    for s in strip.spines.values():
        s.set_color(FRAME); s.set_linewidth(LW_FRAME)
    strip.text(-0.5, 0.5, f"{sv.min():.0f}", ha='left', va='center', fontsize=FS_ANNOT, color=INK)
    strip.text(len(order) - 0.5, 0.5, f"{sv.max():.0f} µg/cm²", ha='right', va='center',
               fontsize=FS_ANNOT, color=INK, fontweight='bold')
    strip.text(nb / 2, -0.25, "NS", ha='center', va='top', fontsize=FS_TICK,
               color=NS_COL, fontweight='bold')
    strip.text((nb + len(order)) / 2, -0.25, "SS", ha='center', va='top',
               fontsize=FS_TICK, color=SS_COL, fontweight='bold')

    # ━━━━━ (B) 로그비 ~ 피지 실측
    frame(b)
    group_scatter(b, md['Sebum'], SL, grp, hollow=LAW0)
    b.set_xlim(-10, 262)
    b.set_xlabel("Sebum (µg/cm²)", color=SEBUM)
    b.set_ylabel(f"log({ital('Staphylococcus')} /\n{ital('Lawsonella')})")   # v8d: 두 줄 — 축 높이보다 길다
    b.set_ylim(-0.4, 10.4)
    head(b, None, None, stat(rho=-0.761, p=0.0001))

    # ━━━━━ (C) 같은 로그비 ~ 수분 — 음성대조
    frame(c)
    group_scatter(c, mr, SL, grp, hollow=LAW0)
    c.set_xlabel("Hydration\n(age-adjusted rank residual)", color=MOIST)   # 두 줄 — 패널 폭을 좁게
    c.tick_params(labelleft=False)
    head(c, None, None, stat(rho=-0.144, p=0.481), sub_col=MUTED)

    # ━━━━━ (D) 속 단위 — 95% CI (검출 기준 통과 4속, 나이 통제, 채움 = BH q < 0.05)
    frame(d, grid='x')
    gen = (GCI[GCI.axis == 'Sebum'].sort_values('rho', ascending=False)['genus'].tolist())
    y0 = np.arange(len(gen))[::-1].astype(float)
    d.axvline(0, color=INK, lw=1.3, zorder=2)
    for yi, g in zip(y0, gen):
        for axis, col, dy in (('Sebum', SEBUM, 0.17), ('Hydration', MOIST, -0.17)):
            r = GCI[(GCI.genus == g) & (GCI.axis == axis)].iloc[0]
            ci_point(d, yi + dy, r.rho, r.lo, r.hi, col, r.q < 0.05, ms=12)
            d.text(0.985, yi + dy, num(r.rho), transform=d.get_yaxis_transform(), ha='right',
                   va='center', fontsize=FS_ANNOT, color=col if r.q < 0.05 else MUTED,
                   fontweight='bold' if r.q < 0.05 else 'normal')
    d.set_yticks(y0); d.set_yticklabels([ital(g) for g in gen])
    d.set_xlim(-1.0, 1.62)   # v8d: 글자를 키워 수치 칸을 넓힘; d.set_ylim(-0.6, len(gen) - 0.4)
    d.set_xticks([-0.8, 0, 0.8])   # v8c: 글자를 키워 −0.4·0.4 를 뺐다 (붙어 보임)
    d.set_xticklabels([num(v, 1) if v else "0" for v in (-0.8, 0, 0.8)])
    d.tick_params(axis='y', length=0)
    d.set_xlabel("partial $\\it{\\rho}$, 95% CI")   # 짧게 — 포스터에서 오른쪽 범례 칩과 겹침

    hs = axis_handles(sebum_first=True)
    fig.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, -0.03), ncol=3,
               frameon=False, fontsize=FS_ANNOT, handletextpad=0.4, columnspacing=1.6,
               borderpad=0)
    group_legend(fig, x=0.985, y=1.0)   # v8b: 빈 원 뜻은 포스터 범례로
    save_fig(fig, name)
    save_panels(fig, {'a': (a, strip), 'b': b, 'c': c, 'd': d}, "F3")
    return fig


if __name__ == "__main__":
    build()
