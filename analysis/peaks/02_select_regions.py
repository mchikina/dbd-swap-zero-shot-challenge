"""Pass 2: pick 10 non-overlapping 2 kb regions, each a strong peak in at least one experiment.
Score per 100 bp bin per experiment = log2((CPM+pseudo)/(median CPM+pseudo)), CPM from replicate-pooled counts;
experiments with < 1M pooled reads are not eligible to nominate a region; bins where a single base holds >30% of
the bin counts are rejected as pile-up artifacts.
Greedy: rank bins by max-over-experiments score, take the top bin, skip anything within 5 kb of
an already chosen region or whose top experiment already contributed a region (so the 10 regions
come from 10 different experiments). chrM and the rDNA repeat (chrXII:451-469 kb) are excluded."""
import numpy as np, pandas as pd, re
d=np.load("bins100.npz"); names=d["names"]; tot=d["total"]; B=d["bins"].astype(np.float64); BX=d["binmax"].astype(np.float64); BIN=int(d["binsize"])
samples=pd.read_csv("/home/mchikina/work/IDR/GSE179430/metadata/samples.tsv",sep="\t").set_index("file").loc[names]
exp=samples["expname"].values
exps=sorted(set(exp), key=str.lower)
# pool replicate counts, then CPM  (138 x NB)
C=np.stack([B[exp==e].sum(0) for e in exps]); CX=np.stack([BX[exp==e].max(0) for e in exps])
depth=np.array([tot[exp==e].sum() for e in exps])
E=C/depth[:,None]*1e6
MINDEPTH=1_000_000
# smooth 3 bins (300 bp) to favour real peaks over single-bin spikes
Es=(E+np.roll(E,1,1)+np.roll(E,-1,1))/3
# chromosome layout (GEO concatenation order)
cs=pd.read_csv("/home/mchikina/work/IDR/genome/sacCer3.chrom.sizes.geo_order",sep="\t",header=None,names=["chr","len"])
cs["start"]=np.concatenate([[0],cs["len"].cumsum().values[:-1]])
binchr=np.empty(E.shape[1],dtype=object); binpos=np.zeros(E.shape[1],dtype=np.int64)
for _,r in cs.iterrows():
    b0=r.start//BIN; b1=(r.start+r.len-1)//BIN
    binchr[b0:b1+1]=r.chr; binpos[b0:b1+1]=np.arange(b0,b1+1)*BIN-r.start
mask=np.ones(E.shape[1],bool)
mask[binchr=="chrM"]=False
mask[(binchr=="chrXII")&(binpos>=451000)&(binpos<=469000)]=False
for _,r in cs.iterrows():  # drop 5 kb at chromosome ends (telomeric repeats)
    mask[(binchr==r.chr)&((binpos<5000)|(binpos>r.len-5000))]=False
pseudo=1.0
med=np.median(Es[:,mask],axis=1,keepdims=True)
S=np.log2((Es+pseudo)/(med+pseudo)); S[:,~mask]=-np.inf
S[depth<MINDEPTH,:]=-np.inf
# reject bins where one base carries >30% of the bin's counts in the top experiment (pile-up artifacts)
spike=CX/np.maximum(C,1)>0.3
top=S.argmax(0); best=S.max(0)
order=np.argsort(-best)
W=2000; MINDIST=5000
chosen=[]; used_exp=set()
for b in order:
    if len(chosen)==10: break
    e=exps[top[b]]
    if e in used_exp or spike[top[b],b]: continue
    c=binchr[b]; center=binpos[b]+BIN//2
    if any(cc==c and abs(center-cen)<MINDIST for cc,cen,*_ in chosen): continue
    start=max(0,center-W//2); start=min(start, int(cs.set_index("chr").loc[c,"len"])-W)
    chosen.append((c,center,start,start+W,e,float(best[b]),float(Es[top[b],b]))); used_exp.add(e)
reg=pd.DataFrame(chosen,columns=["chr","peak_center","start","end","top_experiment","log2_fold_over_median","peak_cpm_per_100bp"])
reg["top_experiment_pooled_depth"]=[int(depth[exps.index(e)]) for e in reg.top_experiment]
reg.insert(0,"region",[f"R{i+1}" for i in range(len(reg))])
# genome-concatenated offsets for window extraction
off=cs.set_index("chr")["start"]
reg["gstart"]=[off[c]+s for c,s in zip(reg.chr,reg.start)]
reg.to_csv("regions.tsv",sep="\t",index=False)
print(reg.to_string())
