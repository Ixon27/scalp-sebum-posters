import numpy as np, math
from xlsx_loader import load
from scipy import stats
rng = np.random.default_rng(20260825)
d = load('../../data/PRJNA1199187/fmicb-2025-1595030-TableS1.xlsx')

def frame(sheet):
    g = d[sheet]; hdr = g[0][1:]
    ids = [r[0] for r in g[1:]]
    M = np.array([[float(x) if x not in ('', None) else 0.0 for x in r[1:1+len(hdr)]] for r in g[1:]])
    return ids, hdr, M

bid, bcols, B = frame('top30.sample.g')
fid, fcols, F = frame('Fungi_top30.sample.g')
meta = {r[0]: r for r in d['metadata'][1:]}
print("세균표", B.shape, "곰팡이표", F.shape)
print("행합 세균 min/max %.4f/%.4f  곰팡이 %.4f/%.4f" % (B.sum(1).min(), B.sum(1).max(), F.sum(1).min(), F.sum(1).max()))

# 두 표의 시료 순서 정렬
order = [fid.index(s) for s in bid]
F = F[order]; fid = [fid[i] for i in order]
assert bid == fid
ids = bid
city = np.array([1 if meta[s][3] == 'Shantou' else 0 for s in ids])
print("시료", len(ids), "| Shantou", int(city.sum()), "Putian", int((1-city).sum()))

def shannon(M):
    P = M / M.sum(1, keepdims=True)
    return np.array([-(p[p > 0] * np.log(p[p > 0])).sum() for p in P])

def get(cols, M, name):
    return M[:, cols.index(name)]

bSh = shannon(B); fSh = shannon(F)
bDom = B.max(1) / B.sum(1)

# 발표된 alpha 와 대조 (21개 시료만 공개돼 있다)
ab = d['alpha_index bacteria']; apub = {r[0]: float(r[6]) for r in ab[1:]}
common = [s for s in ids if s in apub]
if common:
    idx = [ids.index(s) for s in common]
    r, p = stats.spearmanr([apub[s] for s in common], bSh[idx])
    print(f"\n[검증] 발표 Shannon vs 상위30 절단 Shannon (n={len(common)}): rho={r:+.3f}, p={p:.2g}")

# 로그비 — 0 은 검출한계 대체(비영 최소값의 1/2)
def logratio(num, den):
    a = get(bcols, B, num).copy(); b = get(bcols, B, den).copy()
    pc = min(a[a > 0].min(), b[b > 0].min()) / 2
    return np.log((a + pc) / (b + pc)), pc

LR_SL, pc1 = logratio('Staphylococcus', 'Lawsonella')
LR_CL, pc2 = logratio('Cutibacterium', 'Lawsonella')
print(f"\n로그비 pseudocount: Staph/Laws {pc1:.2e}, Cuti/Laws {pc2:.2e}")
for nm in ('Cutibacterium','Staphylococcus','Lawsonella'):
    v = get(bcols, B, nm); print(f"  {nm:16s} 평균 {v.mean()*100:5.1f}%  0인 시료 {int((v==0).sum())}")
Mal = get(fcols, F, 'Malassezia')
print(f"  Malassezia(ITS)  평균 {Mal.mean()*100:5.1f}%  중앙값 {np.median(Mal)*100:5.1f}%  최소 {Mal.min()*100:.1f}%  0인 시료 {int((Mal==0).sum())}")

def perm_spearman(x, y, n=9999, covar=None, label=""):
    if covar is not None:
        def resid(v):
            X = np.c_[np.ones_like(covar, dtype=float), covar.astype(float)]
            return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
        x = resid(stats.rankdata(x)); y = resid(stats.rankdata(y))
        obs = stats.pearsonr(x, y)[0]
        null = np.array([stats.pearsonr(rng.permutation(x), y)[0] for _ in range(n)])
    else:
        obs = stats.spearmanr(x, y)[0]
        null = np.array([stats.spearmanr(rng.permutation(x), y)[0] for _ in range(n)])
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (n + 1)
    print(f"  {label:52s} rho={obs:+.3f}  p={p:.4f}")
    return obs, p

print("\n=== 검정 1. 지질 균형 축 ~ 다양성 축 직교성 (도시 통제) ===")
for nm, lr in (('log(Staph/Lawsonella)', LR_SL), ('log(Cuti/Lawsonella)', LR_CL)):
    perm_spearman(lr, bSh, covar=city, label=f"{nm} ~ 세균 Shannon")
perm_spearman(LR_SL, bDom, covar=city, label="log(Staph/Lawsonella) ~ 우점도")

print("\n=== 검정 2. 세균 지질 축 ~ 곰팡이 (16S↔ITS 짝) ===")
for nm, lr in (('log(Staph/Lawsonella)', LR_SL), ('log(Cuti/Lawsonella)', LR_CL)):
    perm_spearman(lr, Mal, covar=city, label=f"{nm} ~ Malassezia 상대존재량")
    perm_spearman(lr, fSh, covar=city, label=f"{nm} ~ 곰팡이 Shannon")
perm_spearman(bSh, fSh, covar=city, label="세균 Shannon ~ 곰팡이 Shannon")

print("\n=== 참고. 도시 효과 ===")
for nm, v in (('log(Staph/Laws)', LR_SL), ('세균 Shannon', bSh), ('Malassezia', Mal), ('곰팡이 Shannon', fSh)):
    u = stats.mannwhitneyu(v[city == 1], v[city == 0])
    print(f"  {nm:18s} Shantou {np.median(v[city==1]):+.3f} vs Putian {np.median(v[city==0]):+.3f}  p={u.pvalue:.4f}")
