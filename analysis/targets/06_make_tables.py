"""Render replicate and mean-track correlation tables as markdown (correlation_tables.md)."""
import pandas as pd
rep=pd.read_csv("replicate_correlation.tsv",sep="\t"); cor=pd.read_csv("mean_track_correlation.tsv",sep="\t"); man=pd.read_csv("manifest.tsv",sep="\t")
fg=set(man[man["all_peak_r_above_0.2"]].set)
out=[]
out.append("## Replicate correlation\n")
out.append("Pearson r on log2(CPM per 100 bp + 1), chrM excluded. `pairwise` = mean over all replicate pairs; `vs mean of others` = mean over replicates of r(replicate, mean of remaining replicates). `peaks` = union of top-1000 bins of the set's 4 mean tracks.\n")
out.append("| set | role | experiment | n rep | pairwise r genome | pairwise r peaks | vs-mean r genome | vs-mean r peaks |")
out.append("|---|---|---|---|---|---|---|---|")
pw=rep[rep.rep_j!="mean_of_others"].groupby("expname")[["r_genome","r_peaks"]].mean()
vm=rep[rep.rep_j=="mean_of_others"].groupby("expname")[["r_genome","r_peaks"]].mean()
for _,m in man.iterrows():
    tag=" (foreground)" if m.set in fg else ""
    out.append(f"| {m.set}{tag} | {m.role} | {m.expname} | {m.n_rep} | {pw.loc[m.expname,'r_genome']:.3f} | {pw.loc[m.expname,'r_peaks']:.3f} | {vm.loc[m.expname,'r_genome']:.3f} | {vm.loc[m.expname,'r_peaks']:.3f} |")
out.append("\n## Correlation of mean tracks across the 4 experiments of each set\n")
out.append("Pearson r between replicate-mean tracks, log2(CPM per 100 bp + 1). A = wild-type TF A, B = wild-type TF B, A_B_DBD = A body carrying B's DNA-binding domain, B_A_DBD = B body carrying A's DNA-binding domain.\n")
cols=[("A","B"),("A","A_B_DBD"),("A","B_A_DBD"),("B","A_B_DBD"),("B","B_A_DBD"),("A_B_DBD","B_A_DBD")]
for metric,title in [("r_genome","Genome-wide (all 100 bp bins, chrM excluded)"),("r_peaks","At peaks (union of top-1000 bins of the 4 mean tracks)")]:
    out.append(f"### {title}\n")
    out.append("| set | "+" | ".join(f"{a} vs {b}" for a,b in cols)+" | n peak bins |")
    out.append("|---|"+"---|"*(len(cols)+1))
    for s in man.set.unique():
        sub=cor[cor.set==s]; tag=" (foreground)" if s in fg else ""
        vals=[]
        for a,b in cols:
            v=sub[((sub.role_1==a)&(sub.role_2==b))|((sub.role_1==b)&(sub.role_2==a))][metric].iloc[0]; vals.append(f"{v:.3f}")
        out.append(f"| {s}{tag} | "+" | ".join(vals)+f" | {int(sub.n_peak_bins.iloc[0])} |")
    out.append("")
open("correlation_tables.md","w").write("\n".join(out)); print("\n".join(out))
