import numpy as np
from xlsx_loader import load
from scipy import stats
rng=np.random.default_rng(20260825)
d=load('../../data/PRJNA1199187/fmicb-2025-1595030-TableS1.xlsx')
def frame(sh):
    g=d[sh]; hdr=g[0][1:]; ids=[r[0] for r in g[1:]]
    M=np.array([[float(x) if x not in ('',None) else 0.0 for x in r[1:1+len(hdr)]] for r in g[1:]])
    return ids,hdr,M
bid,bcols,B=frame('top30.sample.g'); fid,fcols,F=frame('Fungi_top30.sample.g')
F=F[[fid.index(s) for s in bid]]; ids=bid
meta={r[0]:r for r in d['metadata'][1:]}
city=np.array([1 if meta[s][3]=='Shantou' else 0 for s in ids])
def sh_(M):
    P=M/M.sum(1,keepdims=True)
    return np.array([-(p[p>0]*np.log(p[p>0])).sum() for p in P])
def get(c,M,n): return M[:,c.index(n)]
def pspear(x,y,covar,label,n=9999):
    def resid(v):
        X=np.c_[np.ones_like(covar,dtype=float),covar.astype(float)]
        return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
    xr=resid(stats.rankdata(x)); yr=resid(stats.rankdata(y))
    obs=stats.pearsonr(xr,yr)[0]
    null=np.array([stats.pearsonr(rng.permutation(xr),yr)[0] for _ in range(n)])
    p=(np.sum(np.abs(null)>=abs(obs))+1)/(n+1)
    print(f"  {label:56s} rho={obs:+.3f}  p={p:.4f}  n={len(x)}")
    return obs,p
pc=min(get(bcols,B,'Staphylococcus')[get(bcols,B,'Staphylococcus')>0].min(),
       get(bcols,B,'Lawsonella')[get(bcols,B,'Lawsonella')>0].min())/2
LR_SL=np.log((get(bcols,B,'Staphylococcus')+pc)/(get(bcols,B,'Lawsonella')+pc))
LR_CL=np.log((get(bcols,B,'Cutibacterium')+pc)/(get(bcols,B,'Lawsonella')+pc))
bSh=sh_(B)

print("=== 검정 1a. 발표된 전체깊이 Shannon 으로 재검 (공개된 21개 시료만) ===")
ab=d['alpha_index bacteria']; apub={r[0]:float(r[6]) for r in ab[1:]}
com=[s for s in ids if s in apub]; idx=[ids.index(s) for s in com]
pub=np.array([apub[s] for s in com])
print(f"  (대조) 절단 Shannon vs 발표 Shannon                      rho={stats.spearmanr(bSh[idx],pub)[0]:+.3f}")
for nm,lr in (('log(Staph/Laws)',LR_SL),('log(Cuti/Laws)',LR_CL)):
    pspear(lr[idx],pub,city[idx],f"{nm} ~ 발표 Shannon(전체 깊이)")
    pspear(lr[idx],bSh[idx],city[idx],f"{nm} ~ 절단 Shannon(같은 21개)")

print("\n=== 검정 1b. 3속을 뺀 나머지 28속으로 다양성 재계산 (수학적 얽힘 제거) ===")
keep=[i for i,c in enumerate(bcols) if c not in ('Cutibacterium','Staphylococcus','Lawsonella','Others')]
R=B[:,keep]; ok=R.sum(1)>0
print(f"  나머지 28속 합계: 평균 {R.sum(1).mean()*100:.2f}%  0인 시료 {int((~ok).sum())}개(제외)")
rSh=np.full(len(ids),np.nan); rSh[ok]=sh_(R[ok])
for nm,lr in (('log(Staph/Laws)',LR_SL),('log(Cuti/Laws)',LR_CL)):
    pspear(lr[ok],rSh[ok],city[ok],f"{nm} ~ 잔여 28속 Shannon")

print("\n=== 참고. 우리 코호트와 우점 구조 비교 ===")
print(f"  이 코호트: Cutibacterium {get(bcols,B,'Cutibacterium').mean()*100:.1f}% > Staphylococcus {get(bcols,B,'Staphylococcus').mean()*100:.1f}% > Lawsonella {get(bcols,B,'Lawsonella').mean()*100:.1f}%  (3속 합 {(get(bcols,B,'Cutibacterium')+get(bcols,B,'Staphylococcus')+get(bcols,B,'Lawsonella')).mean()*100:.1f}%)")
print( "  우리 27명 : Staphylococcus 69.5% > Cutibacterium 14.4% > Lawsonella 4.8%  (3속 합 88.7%)")
