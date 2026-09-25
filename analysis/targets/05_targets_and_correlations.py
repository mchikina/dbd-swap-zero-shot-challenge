"""Mean tracks for the 7 paralog pairs with complete DBD-swap sets (A, B, A_B_DBD, B_A_DBD), plus correlations.

Mean track = mean over replicates of per-base CPM (each replicate normalised to its own depth), float32 .npy,
12,157,105 values in GEO order (chrI..chrXVI, chrM).

Correlations use 100 bp bins, log2(CPM per bin + 1), Pearson r, excluding chrM.
Peaks per set = union over the 4 mean tracks of each track's top 1000 bins (excluding chrM, rDNA, 5 kb chromosome ends,
and bins where a single base carries >30% of the bin counts in that experiment).
"""
import numpy as np, pandas as pd, os, itertools
from multiprocessing import Pool
RAW="/home/mchikina/work/IDR/GSE179430/raw"; OUT="mean_tracks"
G=12157105; BIN=100
SETS=[("Dot6","Tod6"),("Ecm22","Upc2"),("Fkh1","Fkh2"),("Gis1","Rph1"),("Pdr8","Yrr1"),("Vhr1","Vhr2"),("Yhp1","Yox1")]
samples=pd.read_csv("/home/mchikina/work/IDR/GSE179430/metadata/samples.tsv",sep="\t")
d=np.load("../peaks/bins100.npz"); bnames=list(d["names"]); btot=d["total"]; B=d["bins"].astype(np.float64); BX=d["binmax"].astype(np.float64)
bidx={n:i for i,n in enumerate(bnames)}
cs=pd.read_csv("/home/mchikina/work/IDR/genome/sacCer3.chrom.sizes.geo_order",sep="\t",header=None,names=["chr","len"])
cs["start"]=np.concatenate([[0],cs["len"].cumsum().values[:-1]])
NB=B.shape[1]; binchr=np.empty(NB,dtype=object); binpos=np.zeros(NB,dtype=np.int64)
for _,r in cs.iterrows():
    b0=r.start//BIN; b1=(r.start+r.len-1)//BIN; binchr[b0:b1+1]=r.chr; binpos[b0:b1+1]=np.arange(b0,b1+1)*BIN-r.start
nuc=binchr!="chrM"
peakable=nuc.copy(); peakable[(binchr=="chrXII")&(binpos>=451000)&(binpos<=469000)]=False
for _,r in cs.iterrows(): peakable[(binchr==r.chr)&((binpos<5000)|(binpos>r.len-5000))]=False

# ---- manifest
rows=[]
for a,b in SETS:
    for role,e in [("A",a),("B",b),("A_B_DBD",f"{a}_{b}_DBD"),("B_A_DBD",f"{b}_{a}_DBD")]:
        sub=samples[samples.expname==e]
        rows.append({"set":f"{a}/{b}","role":role,"expname":e,"n_rep":len(sub),"files":";".join(sub.file),
                     "total_counts":";".join(str(btot[bidx[f]]) for f in sub.file),"mean_track":f"{OUT}/{e}.mean_cpm.npy"})
man=pd.DataFrame(rows); man.to_csv("manifest.tsv",sep="\t",index=False)

# ---- mean tracks (base resolution)
def load(f): return pd.read_csv(f"{RAW}/{f}",header=None,dtype=np.int32).values.ravel()
def mean_track(e):
    fs=samples[samples.expname==e].file.tolist()
    acc=np.zeros(G,dtype=np.float64)
    for f in fs:
        a=load(f); acc+=a/a.sum()*1e6
    acc/=len(fs)
    np.save(f"{OUT}/{e}.mean_cpm.npy",acc.astype(np.float32))
    return e,len(fs),float(acc.sum())
with Pool(8) as p: res=p.map(mean_track,man.expname.tolist())
print("mean tracks written:",len(res))

# ---- binned log CPM per sample and per mean track
L=np.log2(B/btot[:,None]*1e6+1)                          # per-sample log2 CPM per 100bp bin
def r(x,y,m): return float(np.corrcoef(x[m],y[m])[0,1])
mean_bins={}; peaks_set={}
for a,b in SETS:
    exps=man[man.set==f"{a}/{b}"].expname.tolist()
    pk=np.zeros(NB,bool)
    for e in exps:
        idx=[bidx[f] for f in samples[samples.expname==e].file]
        cpm=B[idx]/btot[idx,None]*1e6; m=cpm.mean(0); mean_bins[e]=np.log2(m+1)
        spike=(BX[idx]/np.maximum(B[idx],1)>0.3).any(0)
        ok=peakable&~spike
        top=np.argsort(-np.where(ok,m,-1))[:1000]; pk[top]=True
    peaks_set[f"{a}/{b}"]=pk

# ---- replicate correlations
rr=[]
for _,mrow in man.iterrows():
    e=mrow.expname; pk=peaks_set[mrow.set]
    idx=[bidx[f] for f in samples[samples.expname==e].file]
    for i,j in itertools.combinations(range(len(idx)),2):
        rr.append({"set":mrow.set,"expname":e,"rep_i":i+1,"rep_j":j+1,"r_genome":r(L[idx[i]],L[idx[j]],nuc),"r_peaks":r(L[idx[i]],L[idx[j]],pk)})
    for i in range(len(idx)):
        others=[k for k in idx if k!=idx[i]]
        mo=np.log2((B[others]/btot[others,None]*1e6).mean(0)+1)
        rr.append({"set":mrow.set,"expname":e,"rep_i":i+1,"rep_j":"mean_of_others","r_genome":r(L[idx[i]],mo,nuc),"r_peaks":r(L[idx[i]],mo,pk)})
rep=pd.DataFrame(rr); rep.to_csv("replicate_correlation.tsv",sep="\t",index=False,float_format="%.4f")

# ---- mean-track correlations within each set
cc=[]
for a,b in SETS:
    s=f"{a}/{b}"; sub=man[man.set==s]; pk=peaks_set[s]
    for (r1,e1),(r2,e2) in itertools.combinations(zip(sub.role,sub.expname),2):
        cc.append({"set":s,"role_1":r1,"role_2":r2,"exp_1":e1,"exp_2":e2,"r_genome":r(mean_bins[e1],mean_bins[e2],nuc),
                   "r_peaks":r(mean_bins[e1],mean_bins[e2],pk),"n_peak_bins":int(pk.sum())})
cor=pd.DataFrame(cc); cor.to_csv("mean_track_correlation.tsv",sep="\t",index=False,float_format="%.4f")

# ---- summaries
print("\n== replicate r (pairwise, mean per experiment) ==")
pw=rep[rep.rep_j!="mean_of_others"].groupby(["set","expname"])[["r_genome","r_peaks"]].mean().round(3)
print(pw.to_string())
print("\n== mean-track r within sets ==")
print(cor.pivot_table(index="set",columns=["role_1","role_2"],values="r_genome").round(3).to_string())
print(cor.pivot_table(index="set",columns=["role_1","role_2"],values="r_peaks").round(3).to_string())
