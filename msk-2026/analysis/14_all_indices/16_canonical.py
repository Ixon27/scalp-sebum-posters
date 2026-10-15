"""16  본분석 산출물(01_alpha.tsv·01_logratios.tsv)로 홍반을 다시 잰다 — 초록에 넣을 확정 수치."""
import csv, collections
import numpy as np
from scipy import stats
SEED,NPERM=20260825,9999
rng=np.random.default_rng(SEED)
ROOT = str(__import__("pathlib").Path(__file__).resolve().parents[2])   # = msk-2026/
rep=[]
def W(*a):
    s="".join(str(x) for x in a); rep.append(s); print(s,flush=True)

def load(fn):
    d={}
    with open(f"{ROOT}/output/{fn}") as f:
        r=csv.reader(f,delimiter="\t"); hdr=next(r)[1:]
        for row in r: d[row[0]]=dict(zip(hdr,[float(x) for x in row[1:]]))
    return d,hdr
alpha,acols=load("01_alpha.tsv")
lr,lcols=load("01_logratios.tsv")
md={}
with open(f"{ROOT}/output/01_metadata.tsv") as f:
    r=csv.DictReader(f,delimiter="\t"); idc=r.fieldnames[0]
    for row in r: md[row[idc]]=row
S=[s for s in alpha]
IDX=['Sebum','Moisture','Age','Burnish','Desquamation','Temperature','Redness']
X={k:np.array([float(md[s][k]) for s in S]) for k in IDX}

def pcor(x,y,Z=(),n=NPERM):
    rx,ry=stats.rankdata(x),stats.rankdata(y)
    M=np.c_[np.ones(len(x))] if not Z else np.c_[np.ones(len(x)),*[stats.rankdata(z) for z in Z]]
    f=lambda v: v-M@np.linalg.lstsq(M,v,rcond=None)[0]
    a,b=f(rx),f(ry); o=np.corrcoef(a,b)[0,1]
    null=np.array([np.corrcoef(rng.permutation(a),b)[0,1] for _ in range(n)])
    return o,(np.sum(np.abs(null)>=abs(o))+1)/(n+1)
def bh(ps):
    ps=np.asarray(ps,float); m=len(ps); o=np.argsort(ps); q=np.empty(m); prev=1.0
    for rank,i in enumerate(o[::-1]):
        prev=min(prev,ps[i]*m/(m-rank)); q[i]=prev
    return q

W("="*78); W("16  본분석 산출물 기준 — 초록에 넣을 확정 수치"); W("="*78)
W(f"\n  입력: output/01_alpha.tsv · 01_logratios.tsv (본분석 rarefy 결과) · 시료 {len(S)}개")
W("  ※ 15 번 스크립트는 rarefy 를 새로 뽑아 소수 셋째 자리가 미세하게 달랐다. 보고값은 이 표를 쓴다.")

TAR=[(c,np.array([alpha[s][c] for s in S])) for c in acols]
W("\n### A. 지표 7개 × alpha 5종 (나이 통제, 지표군 안 BH-FDR)")
for tn,y in TAR:
    ps,rs=[],[]
    for k in IDX:
        r,p=pcor(X[k],y,() if k=='Age' else (X['Age'],)); rs.append(r); ps.append(p)
    qs=bh(ps)
    W(f"\n  ── {tn}")
    for k,r,p,q in sorted(zip(IDX,rs,ps,qs),key=lambda t:t[2])[:4]:
        W(f"     {k:<14}rho={r:+.3f}  p={p:.4f}  q={q:.4f}" + ("  **유의**" if q<0.05 else ""))

W("\n### B. 수분 vs 홍반 — 서로 통제")
for tn,y in TAR:
    r1,p1=pcor(X['Moisture'],y,(X['Age'],)); r2,p2=pcor(X['Redness'],y,(X['Age'],))
    r3,p3=pcor(X['Moisture'],y,(X['Age'],X['Redness'])); r4,p4=pcor(X['Redness'],y,(X['Age'],X['Moisture']))
    W(f"  {tn:<10} 수분 {r1:+.3f}(p={p1:.4f}) →홍반통제 {r3:+.3f}(p={p3:.4f}) {'O' if p3<.05 else 'X'}"
      f"   |  홍반 {r2:+.3f}(p={p2:.4f}) →수분통제 {r4:+.3f}(p={p4:.4f}) {'O' if p4<.05 else 'X'}")

W("\n### C. 피지 축은 홍반과 무관한가 (본분석 로그비)")
for c in lcols:
    y=np.array([lr[s][c] for s in S])
    r1,p1=pcor(X['Redness'],y,(X['Age'],))
    r2,p2=pcor(X['Sebum'],y,(X['Age'],X['Redness']))
    W(f"  {c:<38} 홍반 {r1:+.3f}(p={p1:.4f})   피지(홍반통제) {r2:+.3f}(p={p2:.4f})")

W("\n### D. 초록에 쓸 문장용 수치")
sh=np.array([alpha[s]['Shannon'] for s in S]); dm=np.array([alpha[s]['Dominance'] for s in S])
ob=np.array([alpha[s]['Observed'] for s in S])
for nm,y in (("Shannon",sh),("Dominance",dm),("Observed",ob)):
    r,p=pcor(X['Moisture'],y,(X['Age'],)); W(f"  수분 ~ {nm:<10} rho={r:+.2f}  p={p:.3f}")
r,p=pcor(X['Redness'],ob,(X['Age'],X['Moisture'])); W(f"  홍반 ~ Observed (나이+수분 통제) rho={r:+.2f}  p={p:.3f}")
r,p=pcor(X['Redness'],sh,(X['Age'],X['Moisture'])); W(f"  홍반 ~ Shannon  (나이+수분 통제) rho={r:+.2f}  p={p:.3f}")
r,p=pcor(X['Redness'],X['Moisture'],(X['Age'],)); W(f"  홍반 ~ 수분 (나이 통제)          rho={r:+.2f}  p={p:.3f}")
open(f"{ROOT}/output/16_canonical_redness.txt","w").write("\n".join(rep))
print(f"\n저장: {ROOT}/output/16_canonical_redness.txt")
