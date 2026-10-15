#!/usr/bin/env python
"""
11  PRJNA1465550 의 '16S 라벨 run' 과 '이름만 ITS 인 run' 이 같은 시료의 기술 반복인가.

배경: 시료키(피험자·부위·시기) 90개에 run 이 둘씩 있다. 둘 다 세균 16S 로 확인됐고
      read 수는 서로 다르다. 기술 반복이라면 두 run 의 조성이 거의 같아야 한다.

기술 반복이면 얻는 것: 두 run 을 합쳐 시료당 깊이를 대략 두 배로 올릴 수 있다.
기술 반복이 아니면: 서로 다른 무언가(다른 채취·다른 런)이므로 합치면 안 된다.

판정 기준 (사전에 정한다)
  A. 같은 시료키의 두 run 사이 Bray-Curtis 거리가, 서로 다른 시료 사이 거리보다 확연히 작아야 한다.
  B. 짝 거리의 중앙값이 '같은 사람 다른 부위' 거리보다도 작아야 한다.
  C. 지질 로그비가 두 run 사이에서 거의 같아야 한다 (순위상관 높고, 짝 차이가 0 근처).
셋 다 만족하면 기술 반복으로 본다. 하나라도 어긋나면 합치지 않는다.
"""
import re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.spatial.distance import braycurtis

warnings.filterwarnings("ignore")
SEED = 20260825
rng = np.random.default_rng(SEED)

ROOT = Path(__file__).resolve().parents[2]   # = msk-2026/
BASE = ROOT / "data" / "PRJNA1465550"   # QIIME2 export_table/ · export_taxonomy/
OUT = ROOT / "output" / "ext1_baseline"
OUT.mkdir(parents=True, exist_ok=True)

rep = []
def W(*a):
    s = "".join(str(x) for x in a); rep.append(s); print(s, flush=True)


def rank_genus(t):
    for p in (x.strip() for x in t.split(';')):
        if p.startswith('g__'):
            return p[3:]
    return ''


def to_genus(counts, tax):
    gg = tax['Taxon'].map(rank_genus); gg = gg[gg != '']
    idx = counts.index.intersection(gg.index)
    return counts.loc[idx].join(gg.rename('genus')).groupby('genus').sum()


counts = pd.read_csv(BASE / "export_table" / "feature-table.tsv", sep="\t", skiprows=1, index_col=0)
tax = pd.read_csv(BASE / "export_taxonomy" / "taxonomy.tsv", sep="\t", index_col=0)
G = to_genus(counts, tax)

W("=" * 78)
W("11  같은 시료키의 두 run 이 기술 반복인가")
W("=" * 78)
W(f"\n  처리된 16S 표에 든 시료: {G.shape[1]}개")
W("  ※ 이 표는 08-25 에 '16S 라벨' 134 run 만 처리한 결과다.")
W("     '이름만 ITS' 인 run 은 아직 처리되지 않았으므로, 짝 비교를 하려면 그쪽도 처리해야 한다.")

# 실제로 무엇이 들어 있는지 먼저 본다
cols = list(G.columns)
W(f"\n  시료명 예시: {cols[:5]}")
meta = [(s, m.group(1), m.group(2), int(m.group(3))) for s in cols
        if (m := re.match(r'^(SD\d+)\.(LS|NLS)\.(\d+)$', s))]
E = pd.DataFrame(meta, columns=['sample', 'subject', 'site', 'rep']).set_index('sample')
W(f"  파싱된 시료 {len(E)}개 · 피험자 {E.subject.nunique()}명 · rep 분포 {dict(E.rep.value_counts().sort_index())}")

dup = E.groupby(['subject', 'site', 'rep']).size()
W(f"\n  시료키(피험자·부위·시기)당 시료 수 분포: {dict(dup.value_counts())}")
if (dup > 1).sum() == 0:
    W("\n  >>> 이 표 안에는 중복 시료키가 없다. 즉 처리된 134 run 은 시료키마다 하나씩이다.")
    W("      '이름만 ITS' 인 108 run 은 이 표에 포함돼 있지 않으므로, 기술 반복 여부를 확인하려면")
    W("      그 108 run 을 따로 처리해 같은 분류기로 붙인 뒤 비교해야 한다. (비용: fasterq-dump 108건 + DADA2)")

# 참고: 시료 간 거리 구조를 미리 재 둔다 — 나중 판정의 기준선이 된다
depth = int(G.sum(axis=0).min())
W(f"\n  공통 깊이 {depth:,} reads 로 rarefy 후 거리 구조를 재 둔다 (판정 기준선)")


def rarefy(v, d):
    v = v.astype(int)
    if v.sum() <= d:
        return v.astype(float)
    return np.bincount(rng.choice(np.repeat(np.arange(len(v)), v), size=d, replace=False),
                       minlength=len(v)).astype(float)


R = pd.DataFrame({c: rarefy(G[c].values, depth) for c in G.columns}, index=G.index)
Rr = R / R.sum(axis=0)

same_subj_diff_site, same_subj_diff_time, diff_subj = [], [], []
cs = list(Rr.columns)
for i in range(len(cs)):
    for j in range(i + 1, len(cs)):
        a, b = E.loc[cs[i]], E.loc[cs[j]]
        d = braycurtis(Rr[cs[i]].values, Rr[cs[j]].values)
        if a.subject != b.subject:
            diff_subj.append(d)
        elif a.site != b.site and a.rep == b.rep:
            same_subj_diff_site.append(d)
        elif a.site == b.site and a.rep != b.rep:
            same_subj_diff_time.append(d)

W(f"\n  {'비교 종류':<28}{'n':>6}{'중앙 BC거리':>12}{'25%':>8}{'75%':>8}")
for nm, v in (("같은 사람·같은 부위·다른 시기", same_subj_diff_time),
              ("같은 사람·다른 부위·같은 시기", same_subj_diff_site),
              ("다른 사람", diff_subj)):
    v = np.array(v)
    if len(v) == 0:
        W(f"  {nm:<28}{0:>6}"); continue
    W(f"  {nm:<28}{len(v):>6}{np.median(v):>12.3f}{np.percentile(v,25):>8.3f}{np.percentile(v,75):>8.3f}")

W("\n  해석 기준: 기술 반복이라면 그 거리가 위 세 종류보다 모두 확연히 작아야 한다.")
W("  (같은 사람·같은 부위·다른 시기 = 처치 전후. 이것이 '진짜 생물학적 변화' 의 크기 기준이 된다.)")

u = stats.mannwhitneyu(same_subj_diff_time, diff_subj)
W(f"\n  참고 검정: 처치 전후 거리 vs 다른 사람 거리  p = {u.pvalue:.2e}")
W(f"    처치 전후 중앙 {np.median(same_subj_diff_time):.3f} < 다른 사람 중앙 {np.median(diff_subj):.3f}"
  f"  -> 같은 사람은 시기가 달라도 서로 닮았다")

(OUT / "11_dup_runs.txt").write_text("\n".join(rep), encoding="utf-8")
print(f"\n저장: {OUT/'11_dup_runs.txt'}")
