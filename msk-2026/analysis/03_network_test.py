#!/usr/bin/env python
"""
S1-3  공존 네트워크 '붕괴' 주장에 검정을 붙인다.

옛 포스터: NS 77 edges -> SS 34 edges, mean degree 5.1 -> 2.3, density 0.18 -> 0.08,
           modularity 0.47 -> 0.52 를 '붕괴'로 서술. 검정이 없었다.

여기서 두 가지를 확인한다.
 (A) 표본 수 맞춤 귀무분포 — NS 16개에서 11개를 뽑아 다시 만들기를 999회.
     SS(n=11)의 실측값이 이 분포 안에 들어가면 '붕괴'가 아니라 표본 수 차이다.
 (B) 조성 인지 재구성 — CLR + 희소 부분상관(graphical lasso).
     Spearman 상관은 Staphylococcus 가 70%인 조성 데이터에서 가짜 연결을 만든다.

두 검정을 다 통과해야 포스터에 넣는다.
"""
import warnings
from pathlib import Path
import numpy as np, pandas as pd, networkx as nx
from scipy import stats
from sklearn.covariance import GraphicalLassoCV
warnings.filterwarnings("ignore")

RNG = np.random.default_rng(20260825)
NPERM = 999
TOPN = 30            # 옛 포스터와 동일: 상위 30속
RCUT, PCUT = 0.6, 0.05
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

g = pd.read_csv(OUT / "01_genus_counts.tsv", sep="\t", index_col=0)
md = pd.read_csv(OUT / "01_metadata.tsv", sep="\t", index_col=0)
g = g[md.index]
ra = g / g.sum(axis=0) * 100
top = ra.mean(axis=1).sort_values(ascending=False).head(TOPN).index.tolist()

rep = []
W = rep.append


# ------------------------------------------------------------ 네트워크 만들기
def net_spearman(sub_ra):
    """옛 방식: Spearman |r|>=0.6 & p<0.05."""
    X = sub_ra.loc[top].values
    G = nx.Graph(); G.add_nodes_from(top)
    for i in range(len(top)):
        for j in range(i + 1, len(top)):
            a, b = X[i], X[j]
            if np.std(a) == 0 or np.std(b) == 0:
                continue
            r, p = stats.spearmanr(a, b)
            if abs(r) >= RCUT and p < PCUT:
                G.add_edge(top[i], top[j], weight=abs(r), sign=np.sign(r))
    return G


def net_glasso(sub_counts):
    """조성 인지: CLR -> graphical lasso 희소 부분상관."""
    x = sub_counts.loc[top].values.astype(float) + 0.5
    lg = np.log(x)
    Z = (lg - lg.mean(axis=0)).T                      # 시료 x 속
    Z = Z - Z.mean(axis=0)
    sd = Z.std(axis=0); sd[sd == 0] = 1
    Z = Z / sd
    try:
        m = GraphicalLassoCV(alphas=6, max_iter=200, cv=3).fit(Z)
    except Exception:
        return None
    P = m.precision_
    d = np.sqrt(np.diag(P))
    pc = -P / np.outer(d, d)                          # 부분상관
    G = nx.Graph(); G.add_nodes_from(top)
    for i in range(len(top)):
        for j in range(i + 1, len(top)):
            if abs(pc[i, j]) > 1e-6:
                G.add_edge(top[i], top[j], weight=abs(pc[i, j]), sign=np.sign(pc[i, j]))
    return G


def metrics(G):
    n = G.number_of_nodes(); e = G.number_of_edges()
    if e == 0:
        return dict(edges=0, mean_degree=0.0, density=0.0, modularity=0.0, isolates=n)
    comm = nx.community.greedy_modularity_communities(G, weight='weight')
    return dict(edges=e,
                mean_degree=2 * e / n,
                density=nx.density(G),
                modularity=nx.community.modularity(G, comm, weight='weight'),
                isolates=sum(1 for _, d in G.degree() if d == 0))


def null_from_NS(builder, ns_ids, n_take, source):
    """NS에서 n_take개를 뽑아 네트워크를 다시 만드는 것을 NPERM회."""
    rows = []
    for _ in range(NPERM):
        pick = RNG.choice(ns_ids, size=n_take, replace=False)
        G = builder(source(pick))
        if G is not None:
            rows.append(metrics(G))
    return pd.DataFrame(rows)


def report(name, builder, source):
    ns_ids = md.index[md.Group == 'NS'].to_numpy()
    ss_ids = md.index[md.Group == 'SS'].to_numpy()
    G_ns, G_ss = builder(source(ns_ids)), builder(source(ss_ids))
    if G_ns is None or G_ss is None:
        W(f"  {name}: 네트워크 추정 실패"); return
    m_ns, m_ss = metrics(G_ns), metrics(G_ss)

    W(f"\n### {name}")
    W(f"  실측  NS(n=16): edges {m_ns['edges']}, mean degree {m_ns['mean_degree']:.2f}, "
      f"density {m_ns['density']:.3f}, modularity {m_ns['modularity']:.3f}, 고립 {m_ns['isolates']}")
    W(f"  실측  SS(n=11): edges {m_ss['edges']}, mean degree {m_ss['mean_degree']:.2f}, "
      f"density {m_ss['density']:.3f}, modularity {m_ss['modularity']:.3f}, 고립 {m_ss['isolates']}")

    null = null_from_NS(builder, ns_ids, len(ss_ids), source)
    W(f"\n  [표본 수 맞춤 귀무분포] NS에서 11개씩 {len(null)}회 재추출한 네트워크")
    W(f"  {'지표':<13}{'SS 실측':>10}{'NS11 귀무 중앙(2.5–97.5%)':>32}{'p(양측)':>10}   판정")
    verdict = {}
    for k in ['edges', 'mean_degree', 'density', 'modularity']:
        v = null[k].values; obs = m_ss[k]
        lo, hi = np.percentile(v, [2.5, 97.5])
        p_lo = (np.sum(v <= obs) + 1) / (len(v) + 1)
        p_hi = (np.sum(v >= obs) + 1) / (len(v) + 1)
        p = min(1.0, 2 * min(p_lo, p_hi))
        ok = "유의" if p < 0.05 else "설명 안 됨(표본 수로 충분)"
        verdict[k] = (obs, np.median(v), lo, hi, p, p < 0.05)
        W(f"  {k:<13}{obs:>10.3f}{f'{np.median(v):.3f} ({lo:.3f}–{hi:.3f})':>32}{p:>10.4f}   {ok}")
    W(f"\n  -> NS 전체(n=16) 대비 SS의 감소분 중 표본 수만으로 설명되는 부분: "
      f"edges NS16={m_ns['edges']} vs NS11 귀무중앙={verdict['edges'][1]:.1f} vs SS={m_ss['edges']}")
    return G_ns, G_ss, verdict


W("=" * 78)
W("S1-3  공존 네트워크 '붕괴' 검정")
W(f"      상위 {TOPN}속, Spearman |r|>={RCUT} & p<{PCUT} (옛 포스터와 동일 기준)")
W(f"      귀무 재추출 {NPERM}회, 난수고정 20260825")
W("=" * 78)

res_sp = report("(A) Spearman 네트워크 — 옛 포스터 방식", net_spearman, lambda ids: ra[ids])
res_gl = report("(B) 조성 인지 네트워크 — CLR + graphical lasso", net_glasso, lambda ids: g[ids])

# ---- 두 방법이 같은 연결을 보나
if res_sp and res_gl:
    W("\n### 두 방법의 엣지 일치도 (Jaccard)")
    for lab, i in [("NS", 0), ("SS", 1)]:
        A = {frozenset(e) for e in res_sp[i].edges()}
        B = {frozenset(e) for e in res_gl[i].edges()}
        j = len(A & B) / len(A | B) if (A | B) else float('nan')
        W(f"  {lab}: Spearman {len(A)}개 / glasso {len(B)}개 / 공통 {len(A & B)}개  Jaccard = {j:.3f}")

W("\n### 결론")
if res_sp:
    sig = [k for k, v in res_sp[2].items() if v[5]]
    W(f"  Spearman 네트워크에서 표본 수로 설명되지 않는 지표: {sig if sig else '없음'}")
if res_gl:
    sig = [k for k, v in res_gl[2].items() if v[5]]
    W(f"  조성 인지 네트워크에서 표본 수로 설명되지 않는 지표: {sig if sig else '없음'}")

txt = "\n".join(rep)
(OUT / "03_network_test.txt").write_text(txt, encoding="utf-8")
print(txt)
