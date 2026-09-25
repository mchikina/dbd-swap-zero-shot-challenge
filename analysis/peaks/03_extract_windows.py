"""Pass 3: re-read all 348 files and cut out the 10 regions at base resolution."""
import numpy as np, pandas as pd, glob, os
from multiprocessing import Pool
RAW="/home/mchikina/work/IDR/GSE179430/raw"
reg=pd.read_csv("regions.tsv",sep="\t"); W=int(reg.end[0]-reg.start[0])
files=sorted(glob.glob(f"{RAW}/*.out.txt.gz"))
def load(f):
    a=pd.read_csv(f,header=None,dtype=np.int32).values.ravel()
    return os.path.basename(f), a.sum(), np.stack([a[g:g+W] for g in reg.gstart])
with Pool(12) as p: res=p.map(load,files,chunksize=4)
np.savez_compressed("windows.npz",names=np.array([r[0] for r in res]),total=np.array([r[1] for r in res]),
                    win=np.stack([r[2] for r in res]).astype(np.int32),W=W)
print("done", len(res), W)
