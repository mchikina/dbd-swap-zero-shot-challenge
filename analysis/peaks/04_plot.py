"""Plot all experiments across the 10 selected regions.
Fig A: heatmap, rows = 138 experiments (replicate-pooled CPM), columns = 10 regions x 2 kb at 25 bp.
Fig B: line plots per region, every experiment as a thin gray line, top 3 by peak height coloured."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import re

d=np.load("windows.npz"); names=d["names"]; tot=d["total"]; win=d["win"].astype(np.float64); W=int(d["W"])
reg=pd.read_csv("regions.tsv",sep="\t")
samples=pd.read_csv("/home/mchikina/work/IDR/GSE179430/metadata/samples.tsv",sep="\t").set_index("file").loc[names]
exp=samples["expname"].values

# ---- experiment ordering: by paralog pair, then wt / deletion / lactis / DBD
def parse(e):
    m=re.match(r"([A-Za-z0-9]+)(?:_(d([A-Z0-9]+)|lactis|([A-Za-z0-9]+)_DBD))?$",e)
    tf=m.group(1); kind=0; partner=None
    if m.group(2) is None: kind=0
    elif m.group(2).startswith("d"): kind=1; partner=m.group(3).capitalize()
    elif m.group(2)=="lactis": kind=2
    else: kind=3; partner=m.group(4)
    return tf,kind,partner
exps=sorted(set(exp),key=str.lower)
info={e:parse(e) for e in exps}
# pair id from deletion / DBD names
pair={}
for e,(tf,k,p) in info.items():
    if p: 
        key=tuple(sorted([tf.lower(),p.lower()])); pair[tf.lower()]=key; pair[p.lower()]=key
for e,(tf,k,p) in info.items(): pair.setdefault(tf.lower(),(tf.lower(),"zzz"))
exps=sorted(exps,key=lambda e:(pair[info[e][0].lower()],info[e][0].lower(),info[e][1],e))
E=np.stack([win[exp==e].sum(0)/tot[exp==e].sum()*1e6 for e in exps])   # 138 x 10 x W, CPM per base from pooled replicates
nrep={e:int((exp==e).sum()) for e in exps}

# ---- palette (dataviz reference instance): sequential blue ramp, categorical slots 1-3
seq=["#fcfcfb","#cde2fb","#9ec5f4","#6da7ec","#3987e5","#256abf","#184f95","#0d366b"]
cmap=LinearSegmentedColormap.from_list("blue_seq",seq)
cat=["#2a78d6","#eb6834","#1baf7a"]
plt.rcParams.update({"font.size":8,"axes.edgecolor":"#c3c2b7","axes.linewidth":0.6,"xtick.color":"#52514e","ytick.color":"#52514e",
                     "axes.labelcolor":"#0b0b0b","text.color":"#0b0b0b","figure.facecolor":"#fcfcfb","axes.facecolor":"#fcfcfb"})

# ---- Fig A: heatmap
BW=25; nb=W//BW
H=E.reshape(len(exps),len(reg),nb,BW).sum(-1)                  # CPM per 25 bp
vmax=np.percentile(H,99.7); vmin=0.05
fig,axes=plt.subplots(1,len(reg),figsize=(16,20),sharey=True,gridspec_kw={"wspace":0.04})
for i,(ax,r) in enumerate(zip(axes,reg.itertuples())):
    im=ax.imshow(H[:,i,:],aspect="auto",cmap=cmap,norm=LogNorm(vmin=vmin,vmax=vmax),interpolation="nearest",
                 extent=[r.start,r.end,len(exps)-0.5,-0.5])
    ax.set_title(f"{r.region}\n{r.chr}\n{r.start/1000:.1f}-{r.end/1000:.1f} kb\ntop: {r.top_experiment}",fontsize=7)
    ax.set_xticks([r.start+500,r.start+1000,r.start+1500]); ax.set_xticklabels(["-0.5","0","+0.5"],fontsize=6.5)
    ax.tick_params(length=2)
    for s in ax.spines.values(): s.set_visible(False)
axes[0].set_yticks(range(len(exps))); axes[0].set_yticklabels([f"{e} (n={nrep[e]})" for e in exps],fontsize=5.6)
axes[0].tick_params(axis="y",length=0)
# separators between paralog pairs
prev=None
for j,e in enumerate(exps):
    k=pair[info[e][0].lower()]
    if prev is not None and k!=prev:
        for ax in axes: ax.axhline(j-0.5,color="#c3c2b7",lw=0.4)
    prev=k
cb=fig.colorbar(im,ax=axes,fraction=0.012,pad=0.01); cb.set_label("CPM per 25 bp (log scale)",fontsize=8); cb.outline.set_visible(False)
fig.supxlabel("position relative to window centre (kb); window = 2 kb",y=0.005,fontsize=9)
fig.suptitle("GSE179430 ChEC-seq: 138 experiments (replicate-pooled) at 10 peak regions",y=0.995,fontsize=11)
fig.savefig("fig_heatmap_10regions.png",dpi=160,bbox_inches="tight")
fig.savefig("fig_heatmap_10regions.pdf",bbox_inches="tight")
plt.close(fig)

# ---- Fig B: line plots
SM=10                                                          # 10 bp smoothing for lines
def smooth(x): 
    k=np.ones(SM)/SM; return np.convolve(x,k,mode="same")
fig,axes=plt.subplots(5,2,figsize=(13,15)); axes=axes.ravel()
for i,(ax,r) in enumerate(zip(axes,reg.itertuples())):
    x=np.arange(r.start,r.end)
    Y=np.stack([smooth(E[j,i]) for j in range(len(exps))])
    peak=Y.max(1); top3=np.argsort(-peak)[:3]
    for j in range(len(exps)):
        if j in top3: continue
        ax.plot(x,Y[j],color="#c3c2b7",lw=0.5,alpha=0.6)
    for c,j in zip(cat,top3):
        ax.plot(x,Y[j],color=c,lw=1.4,label=f"{exps[j]} (n={nrep[j] if False else nrep[exps[j]]})")
    ax.set_title(f"{r.region}  {r.chr}:{r.start:,}-{r.end:,}",fontsize=9,loc="left")
    ax.legend(frameon=False,fontsize=7,loc="upper right")
    ax.set_ylabel("CPM per bp (10 bp smoothed)",fontsize=7)
    ax.grid(axis="y",color="#eeede9",lw=0.5); ax.set_axisbelow(True)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.ticklabel_format(style="plain",axis="x"); ax.tick_params(labelsize=7)
fig.suptitle("GSE179430 ChEC-seq: all 138 experiments per region (gray), top 3 by peak height coloured",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.98])
fig.savefig("fig_lines_10regions.png",dpi=150); fig.savefig("fig_lines_10regions.pdf")
# table of per-experiment peak height per region (max CPM per 25 bp)
pd.DataFrame(H.max(-1),index=exps,columns=reg.region).round(3).to_csv("peak_height_by_experiment.tsv",sep="\t")
print("plots done")
