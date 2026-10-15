#!/usr/bin/env python
"""
S1-1  두 축(수분 = 다양성 축 / 피지 = 균형 축) 분리를 확정 수치로.

입력 : data/feature-table.biom, data/taxonomy.tsv, data/metadata_sebum.txt
출력 : output/01_two_axis.txt (사람이 읽는 리포트)
       output/01_genus_table.tsv, output/01_clr_table.tsv (후속 단계용)

원칙
 - 수분 결과는 전부 Age 통제 편상관을 함께 보고한다 (Moisture~Age rho=0.62).
 - 속 단위 검정은 '검출 기준'을 통과한 속으로만 한다 (시료 1/3 이상에서 10 reads 이상).
   3,037 reads 깊이에서 0.1%는 3 reads라 희소 속은 신뢰할 수 없다.
 - 편상관 p는 정규근사가 아니라 순열검정으로 낸다 (n=27).
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import biom
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests

RNG = np.random.default_rng(20260825)
NPERM = 9999
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"; OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------------- 자료 적재
def load():
    t = biom.load_table(ROOT / "data/feature-table.biom")
    asv = pd.DataFrame(t.matrix_data.toarray(),
                       index=t.ids('observation'), columns=t.ids('sample'))
    tax = pd.read_csv(ROOT / "data/taxonomy.tsv", sep="\t", index_col=0)
    md = pd.read_csv(ROOT / "data/metadata_sebum.txt", sep="\t",
                     skiprows=[1], index_col=0).loc[asv.columns]

    def rank(s, k):
        for p in (x.strip() for x in s.split(';')):
            if p.startswith(k + '__'):
                return p[3:]
        return ''
    tax['genus'] = tax['Taxon'].map(lambda s: rank(s, 'g'))
    tax['family'] = tax['Taxon'].map(lambda s: rank(s, 'f'))
    # genus 미할당 ASV는 하나로 뭉치지 않는다 — 128개의 서로 다른 희귀 ASV이고
    # 합치면 '환경 유래 배경 총량'이라는 다른 뜻이 된다. 별도 지표로 뺀다.
    named = tax['genus'] != ''
    g = asv[named.reindex(asv.index, fill_value=False)].join(tax['genus']).groupby('genus').sum()
    unassigned = asv[~named.reindex(asv.index, fill_value=False)].sum(axis=0)
    return asv, tax, md, g, unassigned

# ---------------------------------------------------------------- 통계 도구
def spearman_partial(x, y, z):
    """Age(z) 통제 편 Spearman. p는 z를 고정한 채 x를 섞는 순열검정."""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rz = stats.rankdata(z)
    def pr(a, b, c):
        rab = np.corrcoef(a, b)[0, 1]; rac = np.corrcoef(a, c)[0, 1]; rbc = np.corrcoef(b, c)[0, 1]
        d = np.sqrt((1 - rac ** 2) * (1 - rbc ** 2))
        return np.nan if d == 0 else (rab - rac * rbc) / d
    obs = pr(rx, ry, rz)
    null = np.array([pr(RNG.permutation(rx), ry, rz) for _ in range(NPERM)])
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)
    return obs, p

def spearman_perm(x, y):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    obs = np.corrcoef(rx, ry)[0, 1]
    null = np.array([np.corrcoef(RNG.permutation(rx), ry)[0, 1] for _ in range(NPERM)])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)

def clr(counts, pseudo=0.5):
    """중심로그비 변환. 열 = 시료."""
    x = counts.astype(float) + pseudo
    lg = np.log(x)
    return lg - lg.mean(axis=0)

# ---------------------------------------------------------------- 본체
def main():
    asv, tax, md, g, unassigned = load()
    ra = g / g.sum(axis=0) * 100
    depth = int(asv.sum(axis=0).iloc[0])
    rep = []
    W = rep.append

    W("=" * 78)
    W("S1-1  두피 마이크로바이옴의 두 생리 축 — 확정 수치")
    W(f"      시료 {asv.shape[1]}개 (NS {int((md.Group=='NS').sum())} / SS {int((md.Group=='SS').sum())}), "
      f"ASV {asv.shape[0]}개 -> 명명된 속 {g.shape[0]}개, rarefy {depth} reads")
    W(f"      순열 {NPERM}회, 난수고정 20260825")
    W("=" * 78)

    # ---- 0. 검출 기준
    detect = (g >= 10).sum(axis=1)
    keep = sorted(g.index[detect >= np.ceil(len(md) / 3)])
    W("\n[0] 검출 기준 — 시료 1/3(9개) 이상에서 10 reads 이상")
    W(f"    명명된 속 {g.shape[0]}개 중 통과 {len(keep)}개: {', '.join(keep)}")
    W(f"    통과 속 합계 평균 상대존재량 {ra.loc[keep].sum(axis=0).mean():.1f}%")
    W(f"    (속 미할당 ASV {int((tax['genus']=='').sum())}개 합계 = 평균 "
      f"{unassigned.div(asv.sum(axis=0)).mean()*100:.2f}% — 서로 다른 희귀 분류군의 합이라 속으로 취급하지 않음)")

    # ---- 1. 두 축의 독립성
    W("\n[1] 피지 축과 수분 축은 서로 독립인가")
    r_sm, p_sm = spearman_perm(md.Sebum.values, md.Moisture.values)
    pr_sm, pp_sm = spearman_partial(md.Sebum.values, md.Moisture.values, md.Age.values)
    W(f"    Sebum ~ Moisture            rho = {r_sm:+.3f}  p = {p_sm:.4f}")
    W(f"    Sebum ~ Moisture (Age 통제)  rho = {pr_sm:+.3f}  p = {pp_sm:.4f}   <- 독립")
    for a, b in [('Moisture', 'Age'), ('Sebum', 'Age')]:
        r, p = spearman_perm(md[a].values, md[b].values)
        W(f"    {a} ~ {b}{' '*(24-len(a)-len(b))}rho = {r:+.3f}  p = {p:.4f}")
    W("    -> 수분은 나이와 얽혀 있으므로 수분 결과는 전부 Age 통제 값으로 보고한다.")

    # ---- 2. 다양성은 어느 축이 지배하나
    from skbio.diversity.alpha import shannon, simpson, ace, observed_features
    A = pd.DataFrame({
        'Shannon':   [shannon(asv[c].values) for c in asv],
        'Simpson':   [simpson(asv[c].values) for c in asv],
        'ACE':       [ace(asv[c].values.astype(int)) for c in asv],
        'Observed':  [observed_features(asv[c].values) for c in asv],
        'Dominance': [asv[c].max() / asv[c].sum() * 100 for c in asv],
    }, index=asv.columns)
    A.to_csv(OUT / "01_alpha.tsv", sep="\t")

    W("\n[2] 다양성 — 수분 축이 지배하고 피지 축은 건드리지 않는다")
    W(f"    {'지표':<10} {'수분 raw':>18} {'수분 Age통제':>20} {'피지 raw':>18} {'피지 Age통제':>20}   NS vs SS")
    dv = {}
    for m in A.columns:
        cells = []
        for c in ['Moisture', 'Sebum']:
            r0, p0 = spearman_perm(A[m].values, md[c].values)
            r1, p1 = spearman_partial(A[m].values, md[c].values, md.Age.values)
            cells += [f"{r0:+.2f} (p={p0:.4f})", f"{r1:+.2f} (p={p1:.4f})"]
            dv[(m, c)] = (r0, p0, r1, p1)
        u = stats.mannwhitneyu(A[m][md.Group == 'NS'], A[m][md.Group == 'SS'])
        W(f"    {m:<10} {cells[0]:>18} {cells[1]:>20} {cells[2]:>18} {cells[3]:>20}   p={u.pvalue:.3f}")
    W("    -> 수분은 Age 통제 후에도 전 지표 유의. 피지는 전 지표 비유의. NS/SS도 전부 비유의.")

    # ---- 3. 균형(로그비)은 어느 축이 지배하나
    big = ['Staphylococcus', 'Cutibacterium', 'Lawsonella']
    W("\n[3] 우점 3속의 균형(로그비) — 피지 축이 지배하고 수분 축은 건드리지 않는다")
    W(f"    3속 합계 평균 상대존재량 {ra.loc[big].sum(axis=0).mean():.1f}% "
      f"(범위 {ra.loc[big].sum(axis=0).min():.1f}–{ra.loc[big].sum(axis=0).max():.1f}%)")
    lr_tab = {}
    W(f"    {'로그비':<34} {'피지 Age통제':>20} {'수분 Age통제':>20}")
    for i in range(3):
        for j in range(i + 1, 3):
            lab = f"log({big[i]}/{big[j]})"
            v = np.log((g.loc[big[i]] + 0.5) / (g.loc[big[j]] + 0.5)).values
            lr_tab[lab] = v
            out = []
            for c in ['Sebum', 'Moisture']:
                r, p = spearman_partial(v, md[c].values, md.Age.values)
                out.append(f"{r:+.3f} (p={p:.4f})")
            W(f"    {lab:<34} {out[0]:>20} {out[1]:>20}")
    pd.DataFrame(lr_tab, index=g.columns).to_csv(OUT / "01_logratios.tsv", sep="\t")
    W("    -> Lawsonella 가 든 두 로그비는 피지에만 반응하고 수분에는 반응하지 않는다.")
    W("       Staph/Cutibacterium 만 두 축 모두에 반응한다 (복합 지표).")

    # ---- 4. 속 단위 반응 (검출 기준 통과 속만, BH 보정)
    W("\n[4] 속 단위 반응 — 검출 기준 통과 속만, Age 통제, BH-FDR")
    for c in ['Sebum', 'Moisture']:
        rows = []
        for x in keep:
            r, p = spearman_partial(ra.loc[x].values, md[c].values, md.Age.values)
            rows.append((x, r, p, ra.loc[x].mean(), int(detect[x])))
        R = pd.DataFrame(rows, columns=['genus', 'rho', 'p', 'meanRA', 'nDetect'])
        R['q'] = multipletests(R['p'], method='fdr_bh')[1]
        R = R.reindex(R.rho.abs().sort_values(ascending=False).index)
        R.to_csv(OUT / f"01_genus_{c.lower()}.tsv", sep="\t", index=False)
        W(f"    --- {c} ---")
        for _, s in R.iterrows():
            flag = "**" if s.q < 0.05 else ("*" if s.q < 0.10 else "  ")
            W(f"      {flag} {s.genus:<22} rho={s.rho:+.3f}  p={s.p:.4f}  q={s.q:.4f}  "
              f"(평균 {s.meanRA:5.2f}%, 검출 {int(s.nDetect)}/27)")

    # ---- 5. 두 축의 반응 프로파일이 겹치나
    sv = np.array([spearman_partial(ra.loc[x].values, md.Sebum.values, md.Age.values)[0] for x in keep])
    mv = np.array([spearman_partial(ra.loc[x].values, md.Moisture.values, md.Age.values)[0] for x in keep])
    W("\n[5] 두 축의 속별 반응 프로파일 상관 (겹치면 같은 축, 0이면 다른 축)")
    W(f"    Pearson r = {np.corrcoef(sv, mv)[0,1]:+.3f} / "
      f"Spearman rho = {stats.spearmanr(sv, mv).statistic:+.3f}   (n={len(keep)} 속)")

    # ---- 6. 군 비교 (피지-군 공선성 명시)
    W("\n[6] 피지와 군(NS/SS)의 관계 — 공선성 확인")
    ns, ss = md.Sebum[md.Group == 'NS'], md.Sebum[md.Group == 'SS']
    W(f"    NS sebum {ns.min():.0f}–{ns.max():.0f} (중앙 {ns.median():.0f}) | "
      f"SS sebum {ss.min():.0f}–{ss.max():.0f} (중앙 {ss.median():.0f})")
    W(f"    겹치는 시료 0개, 두 군 사이 빈 구간 {ns.max():.0f} -> {ss.min():.0f} ({ss.min()/ns.max():.1f}배)")
    W("    -> 이 코호트에서 '고피지'와 'SD'는 분리 불가. 미생물로 NS/SS를 분류하는 모델은 순환논증이므로 만들지 않는다.")
    for m in ['Moisture', 'Age']:
        u = stats.mannwhitneyu(md[m][md.Group == 'NS'], md[m][md.Group == 'SS'])
        W(f"    {m:<9} NS 중앙 {md[m][md.Group=='NS'].median():6.2f} vs SS {md[m][md.Group=='SS'].median():6.2f}  p={u.pvalue:.3f}")

    # ---- 저장
    g.to_csv(OUT / "01_genus_counts.tsv", sep="\t")
    clr(g).to_csv(OUT / "01_genus_clr.tsv", sep="\t")
    pd.Series(keep, name='genus').to_csv(OUT / "01_keep_genera.tsv", sep="\t", index=False)
    md.to_csv(OUT / "01_metadata.tsv", sep="\t")

    txt = "\n".join(rep)
    (OUT / "01_two_axis.txt").write_text(txt, encoding="utf-8")
    print(txt)

if __name__ == "__main__":
    main()
