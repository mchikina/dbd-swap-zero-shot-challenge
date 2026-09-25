"""Pass 1: read all 348 GSE179430 per-base coverage files, sum into 100 bp bins."""
import numpy as np, pandas as pd, os, glob
from multiprocessing import Pool
RAW="/home/mchikina/work/IDR/GSE179430/raw"
BIN=100
G=12157105
NB=(G+BIN-1)//BIN
files=sorted(glob.glob(f"{RAW}/*.out.txt.gz"))
def load(f):
    a=pd.read_csv(f,header=None,dtype=np.int32).values.ravel()
    assert a.shape[0]==G, f
    pad=np.zeros(NB*BIN,dtype=np.int64); pad[:G]=a
    r=pad.reshape(NB,BIN)
    return os.path.basename(f), a.sum(), a.max(), r.sum(1).astype(np.int32), r.max(1).astype(np.int32)
with Pool(12) as p:
    res=p.map(load,files,chunksize=4)
names=[r[0] for r in res]
tot=np.array([r[1] for r in res]); mx=np.array([r[2] for r in res])
M=np.stack([r[3] for r in res]); MX=np.stack([r[4] for r in res])
np.savez_compressed("bins100.npz",names=np.array(names),total=tot,maxval=mx,bins=M,binmax=MX,binsize=BIN)
pd.DataFrame({"file":names,"total_counts":tot,"max_per_base":mx}).to_csv("sample_totals.tsv",sep="\t",index=False)
print("done",M.shape, "total counts min/median/max", tot.min(), np.median(tot), tot.max(), "max per-base", mx.max())
