"""
공통 스타일 v7 (2026-10-13 축 색 명도 맞춤·배치 정리; v6 KSABC 팔레트; v4 2026-10-11) — 2026-06 KSABC 포스터(PAM-26) figure 를 기준으로 다시 맞췄다.

━━━ 기준을 어디서 잡았나 ━━━
KSABC 포스터(80 × 120 cm) 파일 `2026KSABC_최종.pptx` 에서 그림마다 배치 너비와 PDF 글자 크기를 재,
포스터에 실제로 찍힌 크기로 환산했다. 그림 안 글씨가 제일 잘 읽힌 F1a·F3b·F5b 의 값이 아래 목표다.
    제목 20–27 pt · 축 제목 18–25 pt · 눈금 15–20 pt · 통계·값 15–21 pt
이 세트는 **포스터 한 단 너비 38.5 cm 에 1:1 로** 그린다. 그 너비로 넣으면 글자가 위 목표 크기로 찍힌다.
(KSABC 는 그림마다 확대율이 1.3–2.3 배로 달라 같은 포스터 안에서 글씨가 11–27 pt 로 들쭉날쭉했다. 이번엔 한 값으로 고정한다.)

━━━ KSABC 에서 그대로 가져온 것 ━━━
1. 색 — 정상 NS = #002060(남색) · 지루성 SS = #C0000D(빨강). 군은 색으로 구분하고 점은 채운 원 + 검은 테두리.
2. 틀 — 네 변을 다 닫은 테두리(#333333) + 연한 회색 주 격자(#EBEBEB). ggplot theme_bw 와 같은 모양.
3. 통계 — "Spearman r = 0.74, p < 0.001" 형식. 통계 기호(ρ·p·R)는 이탤릭, 쉼표로 잇는다, 양수에 + 를 붙이지 않는다.
   관계 하나짜리 패널은 제목 바로 아래 한 줄(부제), 군 비교는 패널 안 위쪽 가운데.
4. 값 라벨 — 마크 옆에 마크와 같은 색, 굵게.
5. 글꼴 — Helvetica 계열 산세리프(여기선 Arial). 속명 이탤릭.

━━━ 이 포스터만의 규칙 (KSABC 에 없던 것) ━━━
6. 축 색 — 수분 = 파랑 #4A8FCF · 피지 = 산호 #D9604F. 군 색(남색·빨강)과 같은 계열의 밝은 톤, 둘의 명도를 맞춤 (v7). 제3 변수(나이·홍반·머신러닝) = 진회색 #404040.
   군 색(남색·빨강)과 겹치지 않게 축에는 다른 색을 쓴다.
7. 비유의 통계는 회색(#7F7F7F). KSABC 의 envfit 그림이 비유의 변수를 회색으로 뺀 방식과 같다. "n.s." 는 붙이지 않는다.
8. 해석 문장은 그림에 넣지 않는다. 그림은 짧은 중립 제목 + 통계값만 갖고, 뜻은 포스터 캡션이 맡는다.
9. 신뢰구간은 가로선, 점 = 추정값. 채운 점 = p < 0.05, 속 빈 점 = 비유의.
   v5 부터 막대·타일(F2D)도 같은 규칙 — 채움 = 유의, 속 빔(흰 면 + 색 테두리) = 비유의.
10. 3속 조성(F3A)은 v4 색(푸른 회색·연어색·진홍) — 사용자 지시로 되돌림.
11. (v8, 2026-10-14) 그림 안에 패널 제목·패널 문자·기호 설명 부제를 넣지 않는다 (KSABC 처럼).
    패널 문자는 포스터에서 텍스트 상자로 붙이고, 기호 뜻(채움 = p < 0.05 등)은 포스터 레전드가 맡는다.
    그림 안 글자는 축 라벨·눈금·통계값·색 범례만. 속명은 줄이지 않는다(Staph·Laws·Cuti 금지).
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

# ━━━ 판 이름 — 다시 그리면 DATE 를 그날로, VER 를 1 올린다 (볼트 규칙 YYMMDD_이름_vN) ━━━
DATE, VER = "261014", 8

# ━━━ 색 — 7장 전체가 이 표만 쓴다. 그림 스크립트에 색 코드를 직접 쓰지 않는다 (v6, 2026-10-12) ━━━
# 파랑·빨강이 주색이다 (사용자 지시 2026-10-13 "파랑, 붉은색이 메인"). 모든 색을 이 두 계열 + 회색에서 뽑는다.
#   군 = 가장 진한 남색 #002060 · 빨강 #C0000D (KSABC 군 색)
#   축 = 같은 계열의 밝고 부드러운 톤 — 수분 파랑 #4A8FCF · 피지 산호 #D9604F.
#        OKLCH 명도 0.63 · 0.64 로 맞췄고, 군 색(0.27 · 0.51)보다 밝아 한 그림에 같이 나와도 갈린다.
#        (v6 하늘색·주황은 명도 0.54 · 0.74 로 어긋났고, v7 초안의 주황은 '붉은 계열'이 아니었다.)
#   제3 변수(나이·홍반·머신러닝) = 진회색 (KSABC F1b envfit 화살표처럼) / 속 조성 = v4 색(푸른 회색·연어색·진홍)
# 회색 글자·선은 "비유의·부수 정보" 뜻으로만 쓴다.
NS_COL = '#002060'   # 정상 (KSABC 그대로)
SS_COL = '#C0000D'   # 지루성 (KSABC 그대로)
SEBUM = '#D9604F'    # 피지 축 — 산호(빨강 계열의 밝은 톤)
MOIST = '#4A8FCF'    # 수분 축 — 파랑
OTHER = '#404040'    # 제3 변수 — 나이(F4)·홍반(F2D)·머신러닝(옛 F5, 삭제). 진회색 (KSABC envfit 화살표)
AGE = OTHER
INK = '#000000'
FRAME = '#333333'    # 틀·박스 선
EDGE = '#1A1A1A'     # 점 테두리
DARK = '#555555'     # 범례의 '모양만' 보이는 견본
MUTED = '#7F7F7F'    # 비유의·부수 정보 (글자·선)
RULE = '#BDBDBD'     # 구분선
PALE = '#E6E6E6'     # 귀무 띠·면
BAND = '#F4F4F4'     # 줄무늬 배경
GRID = '#EBEBEB'
NS_T, SS_T = '#C9D1E3', '#F2C4C6'          # 군 연한 톤 (면·밀도)
SEBUM_T, MOIST_T = '#F6D9D4', '#D3E4F4'    # 축 연한 톤 (띠·신뢰대)
SEBUM_RAMP = ['#FDF2F0', '#EEA99F', SEBUM]  # 피지 실측 띠 (F3A)
GROUP_COL = {'NS': NS_COL, 'SS': SS_COL}

# 3속 조성 막대 — v4 색으로 되돌림 (2026-10-13 사용자 지시 "원래 붉은 것·파란 것").
# Staphylococcus 푸른 회색 · Cutibacterium 연어색 · Lawsonella 진홍. 지루성에서 늘어나는 Lawsonella 가 붉은 계열.
# Cutibacterium 은 v4 #E9A383 에서 한 단계 진하게(OKLCH 명도 0.78 → 0.72) — Staphylococcus 와 명도가 같아
# 적록색맹에서 색조로만 갈리던 것(거리 10)을 명도 차로 벌렸다(13~15). (2026-10-13 사용자 승인)
# (v5 보라 · v6 회색 · v7 초안 초록을 거쳤다.) '기타'는 연회색.
GEN_COL = {'Staphylococcus': '#A9B9CC', 'Cutibacterium': '#DA9167', 'Lawsonella': '#9E1B22',
           'other genera': PALE}

# ━━━ 글자 크기 — 포스터 1:1 (pt) ━━━
FS_TITLE, FS_LAB, FS_TICK, FS_STAT, FS_ANNOT, FS_LETTER = 20, 18, 15, 16, 15, 24
# v8c: 포스터에서 더 많이 줄여 놓는 그림(F2·F3)은 import 전에 FS_SCALE 을 정해 글자를 키운다
_K = float(os.environ.get('FS_SCALE', 1))
FS_TITLE, FS_LAB, FS_TICK, FS_STAT, FS_ANNOT, FS_LETTER = (round(v * _K, 1) for v in
    (FS_TITLE, FS_LAB, FS_TICK, FS_STAT, FS_ANNOT, FS_LETTER))

# ━━━ 선·마크 ━━━
LW_FRAME, LW_GRID, LW_DATA, LW_CI = 1.3, 1.0, 2.6, 2.4
S_PT = 170           # 산점도 점 면적 (지름 약 13 pt)
COL_W = 38.5 / 2.54  # 포스터 한 단 너비 (in)
CM = 1 / 2.54

MINUS = '−'
HALO = [pe.withStroke(linewidth=4, foreground='white')]

mpl.rcParams.update({
    'font.family': 'Arial',
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'font.size': FS_TICK,
    'axes.titlesize': FS_TITLE, 'axes.titleweight': 'bold',
    'axes.labelsize': FS_LAB, 'axes.labelweight': 'bold', 'axes.labelcolor': INK,
    'axes.linewidth': LW_FRAME, 'axes.edgecolor': FRAME,
    'xtick.labelsize': FS_TICK, 'ytick.labelsize': FS_TICK,
    'xtick.major.width': LW_FRAME, 'ytick.major.width': LW_FRAME,
    'xtick.major.size': 6, 'ytick.major.size': 6,
    'xtick.color': FRAME, 'ytick.color': FRAME, 'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
    'legend.fontsize': FS_ANNOT, 'legend.frameon': False,
    'figure.dpi': 72, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[1])   # = msk-2026/
FIGDIR = f"{ROOT}/output/figures"
os.makedirs(f"{FIGDIR}/panels", exist_ok=True)


# ━━━━━━━━━━ 수치 표기 (형식은 여기 한 곳에서만 정한다) ━━━━━━━━━━
def num(v, d=None):
    """상관·계수: 소수 2자리, |v| < 0.1 이면 3자리. 음수만 부호(진짜 마이너스)를 붙인다."""
    if d is None:
        d = 3 if abs(v) < 0.1 else 2
    s = f"{abs(v):.{d}f}"
    return (MINUS + s) if v < 0 and float(s) != 0 else s


def fmt_p(p):
    return "$\\it{p}$ < 0.001" if p < 0.001 else f"$\\it{{p}}$ = {p:.3f}"


def stat(rho=None, p=None, r2=None, kind="partial"):
    """부제 한 줄: 'partial ρ = −0.76, p < 0.001' (KSABC 'Spearman r = 0.74, p < 0.001' 형식)."""
    bits = []
    if rho is not None:
        bits.append((f"{kind} " if kind else "") + f"$\\it{{\\rho}}$ = {num(rho)}")
    if r2 is not None:
        bits.append(f"$\\it{{R}}$² = {r2:.3f}")
    if p is not None:
        bits.append(fmt_p(p))
    return ", ".join(bits)


def ci_text(r, lo, hi):
    return f"{num(r)}  [{num(lo, 2)}, {num(hi, 2)}]"


def sig_col(p, col=INK):
    return col if p < 0.05 else MUTED


def ital(name):
    return f"$\\it{{{name}}}$"


RATIO = f"log({ital('Staphylococcus')} / {ital('Lawsonella')})"
RATIO_S = RATIO   # v8: 속명 줄임 금지
LIPID_S = f"log({ital('Lawsonella')} / {ital('Staphylococcus')})"


# ━━━━━━━━━━ 틀 · 머리 ━━━━━━━━━━
def frame(ax, grid='both'):
    """KSABC theme_bw — 네 변 닫힌 틀 + 연한 주 격자."""
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(FRAME); s.set_linewidth(LW_FRAME)
    ax.tick_params(colors=FRAME, labelcolor=INK, labelsize=FS_TICK)
    if grid:
        ax.set_axisbelow(True)
        ax.grid(axis=grid, color=GRID, lw=LW_GRID, zorder=0)


def head(ax, letter=None, title=None, sub=None, sub_col=INK, lx=-0.12, reserve=True):
    """패널 머리 — 문자(왼쪽 위) · 굵은 제목 · 그 아래 통계 부제. 전부 데이터 영역 밖.
    lx 는 패널 문자의 가로 위치(축 너비 비율). y 눈금 글자 폭에 맞춰 패널마다 준다.
    reserve=True 면 부제가 없어도 그 자리를 비워 둔다 — 같은 행의 제목·문자 높이가 맞는다."""
    title = letter = None   # v8: 제목·패널 문자는 그리지 않는다 (규칙 11)
    reserve = False
    y = 6
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords='axes fraction', xytext=(0, y),
                    textcoords='offset points', ha='left', va='bottom', fontsize=FS_STAT,
                    color=sub_col, annotation_clip=False)
    if sub or reserve:
        y += FS_STAT * 1.45
    if title:
        ax.annotate(title, xy=(0, 1), xycoords='axes fraction', xytext=(0, y),
                    textcoords='offset points', ha='left', va='bottom', fontsize=FS_TITLE,
                    fontweight='bold', color=INK, annotation_clip=False)
    if letter:
        ax.annotate(letter.upper(), xy=(lx, 1), xycoords='axes fraction', xytext=(0, y),
                    textcoords='offset points', ha='left', va='bottom', fontsize=FS_LETTER,
                    fontweight='bold', color=INK, annotation_clip=False)


def inside_top(ax, text, col=INK, x=0.5, ha='center'):
    """군 비교 패널의 p — KSABC F1a 처럼 패널 안 위쪽 가운데."""
    ax.text(x, 0.955, text, transform=ax.transAxes, ha=ha, va='top', fontsize=FS_STAT,
            color=col, zorder=8)


# ━━━━━━━━━━ 마크 ━━━━━━━━━━
def group_scatter(ax, x, y, grp, s=S_PT, hollow=None, zorder=4):
    """군 = 색 (KSABC). hollow 가 True 인 점은 속을 비운다(예: Lawsonella 미검출)."""
    x, y, grp = np.asarray(x, float), np.asarray(y, float), np.asarray(grp)
    hollow = np.zeros(len(x), bool) if hollow is None else np.asarray(hollow, bool)
    for g in ('NS', 'SS'):
        for h in (False, True):
            m = (grp == g) & (hollow == h)
            if m.any():
                ax.scatter(x[m], y[m], s=s, marker='o',
                           facecolor='white' if h else GROUP_COL[g],
                           edgecolor=GROUP_COL[g] if h else EDGE,
                           lw=2.0 if h else 0.9, zorder=zorder, alpha=1.0)


def fit_band(ax, x, y, color, tint, seed=20260825, nboot=2000):
    """선형 적합 + 부트스트랩 95% 신뢰대 (KSABC F5b 처럼 연한 띠 위에 진한 선)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    xs = np.linspace(x.min(), x.max(), 120)
    rng = np.random.default_rng(seed)
    boot = np.full((nboot, xs.size), np.nan)
    for i in range(nboot):
        idx = rng.integers(0, x.size, x.size)
        if np.ptp(x[idx]) > 0:
            boot[i] = np.polyval(np.polyfit(x[idx], y[idx], 1), xs)
    lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    ax.fill_between(xs, lo, hi, color=tint, alpha=0.75, lw=0, zorder=1)
    ax.plot(xs, np.polyval(np.polyfit(x, y, 1), xs), color=color, lw=LW_DATA + 0.6, zorder=2)


def ci_point(ax, y, r, lo, hi, color, filled, marker='o', ms=13, lw=LW_CI, z=4):
    ax.plot([lo, hi], [y, y], color=color, lw=lw, solid_capstyle='butt', zorder=z)
    ax.plot([r], [y], marker=marker, ms=ms, mfc=color if filled else 'white', mec=color,
            mew=2.2, zorder=z + 1, ls='')


def group_legend(fig, x=0.99, y=1.0, extra=None, ncol=None):
    h = [Line2D([], [], marker='o', ls='', ms=13, mfc=NS_COL, mec=EDGE, mew=0.9,
                label="NS"),
         Line2D([], [], marker='o', ls='', ms=13, mfc=SS_COL, mec=EDGE, mew=0.9,
                label="SS")]
    if extra:
        h += extra
    fig.legend(handles=h, loc='upper right', bbox_to_anchor=(x, y), ncol=ncol or len(h),
               frameon=False, fontsize=FS_ANNOT, handletextpad=0.3, columnspacing=1.6,
               borderpad=0)


def axis_handles(filled_note=None, sebum_first=False):
    hs = [Line2D([], [], marker='o', ls='-', lw=LW_CI, ms=12, color=MOIST, mfc=MOIST,
                 label='Hydration'),
          Line2D([], [], marker='o', ls='-', lw=LW_CI, ms=12, color=SEBUM, mfc=SEBUM,
                 label='Sebum')]
    if sebum_first:
        hs = hs[::-1]
    return hs   # v8: 채움/속 빔 뜻은 포스터 레전드에


# ━━━━━━━━━━ 저장 ━━━━━━━━━━
def fname(name):
    return f"{DATE}_{name}_v{VER}"


def save_fig(fig, name):
    base = f"{FIGDIR}/{fname(name)}"
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=300)
    print(f"저장: output/figures/{fname(name)}.pdf / .png")


def save_panels(fig, panels, prefix, pad=0.05):   # v8: 0.25 → 0.05 (포스터에서 패널 간격을 좁게)
    """패널을 하나씩 잘라 따로 저장 — 포스터에서 자유롭게 배치하려고. 머리(문자·제목·부제)까지 포함.
    값이 축 묶음(tuple)이면 함께 저장한다(예: 본 패널 + 주변 밀도). 저장하는 동안 다른 축과
    figure 범례는 숨긴다 — 패널 문자가 이웃 패널 쪽으로 뻗어 있어 잘라낸 범위에 이웃이 걸리기 때문."""
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for tag, axs in panels.items():
        axs = axs if isinstance(axs, (tuple, list)) else (axs,)
        hidden = [o for o in fig.axes if o not in axs and o.get_visible()] + \
                 [lg for lg in fig.legends if lg.get_visible()]
        for o in hidden:
            o.set_visible(False)
        # v8c: 축 폭보다 긴 축 이름은 get_tightbbox 가 잘라 버린다 → 축 이름 범위를 따로 합친다
        ext = [ax.get_tightbbox(r) for ax in axs] + \
              [lab.get_window_extent(r) for ax in axs for lab in (ax.xaxis.label, ax.yaxis.label)
               if lab.get_visible() and lab.get_text()]
        bb = Bbox.union(ext).transformed(fig.dpi_scale_trans.inverted())
        bb = Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
        for ext in ('pdf', 'png'):
            fig.savefig(f"{FIGDIR}/panels/{DATE}_{prefix}{tag}_v{VER}.{ext}", bbox_inches=bb,
                        dpi=300 if ext == 'png' else None)
        for o in hidden:
            o.set_visible(True)
    print(f"  패널 개별 저장: output/figures/panels/{DATE}_{prefix}*_v{VER} ({len(panels)}개)")


# ━━━━━━━━━━ 데이터 (전 figure 공통 입구) ━━━━━━━━━━
def load():
    import pandas as pd
    md = pd.read_csv(f"{ROOT}/output/01_metadata.tsv", sep="\t", index_col=0)
    alpha = pd.read_csv(f"{ROOT}/output/01_alpha.tsv", sep="\t", index_col=0)
    lr = pd.read_csv(f"{ROOT}/output/01_logratios.tsv", sep="\t", index_col=0)
    counts = pd.read_csv(f"{ROOT}/output/01_genus_counts.tsv", sep="\t", index_col=0)
    S = list(md.index)
    return md.loc[S], alpha.loc[S], lr.loc[S], counts[S]


def table(name):
    import pandas as pd
    return pd.read_csv(f"{ROOT}/output/{name}", sep="\t")


def resid(v, *covs):
    """순위 잔차 — 나이(등) 통제 편상관을 그림 좌표로 옮길 때 쓴다."""
    from scipy import stats
    v = np.asarray(v, float)
    Z = np.column_stack([np.ones(v.size)] + [stats.rankdata(c) for c in covs])
    r = stats.rankdata(v)
    return r - Z @ np.linalg.lstsq(Z, r, rcond=None)[0]
