## Replicate correlation

Pearson r on log2(CPM per 100 bp + 1), chrM excluded. `pairwise` = mean over all replicate pairs; `vs mean of others` = mean over replicates of r(replicate, mean of remaining replicates). `peaks` = union of top-1000 bins of the set's 4 mean tracks.

| set | role | experiment | n rep | pairwise r genome | pairwise r peaks | vs-mean r genome | vs-mean r peaks |
|---|---|---|---|---|---|---|---|
| Dot6/Tod6 (foreground) | A | Dot6 | 4 | 0.827 | 0.916 | 0.871 | 0.943 |
| Dot6/Tod6 (foreground) | B | Tod6 | 2 | 0.864 | 0.945 | 0.864 | 0.945 |
| Dot6/Tod6 (foreground) | A_B_DBD | Dot6_Tod6_DBD | 2 | 0.936 | 0.963 | 0.936 | 0.963 |
| Dot6/Tod6 (foreground) | B_A_DBD | Tod6_Dot6_DBD | 3 | 0.495 | 0.581 | 0.567 | 0.622 |
| Ecm22/Upc2 | A | Ecm22 | 2 | 0.900 | 0.950 | 0.900 | 0.950 |
| Ecm22/Upc2 | B | Upc2 | 2 | 0.918 | 0.983 | 0.918 | 0.983 |
| Ecm22/Upc2 | A_B_DBD | Ecm22_Upc2_DBD | 2 | 0.636 | 0.857 | 0.636 | 0.857 |
| Ecm22/Upc2 | B_A_DBD | Upc2_Ecm22_DBD | 4 | 0.808 | 0.868 | 0.862 | 0.909 |
| Fkh1/Fkh2 | A | Fkh1 | 2 | 0.928 | 0.983 | 0.928 | 0.983 |
| Fkh1/Fkh2 | B | Fkh2 | 2 | 0.933 | 0.959 | 0.933 | 0.959 |
| Fkh1/Fkh2 | A_B_DBD | Fkh1_Fkh2_DBD | 4 | 0.939 | 0.955 | 0.958 | 0.970 |
| Fkh1/Fkh2 | B_A_DBD | Fkh2_Fkh1_DBD | 3 | 0.920 | 0.958 | 0.939 | 0.970 |
| Gis1/Rph1 | A | Gis1 | 2 | 0.947 | 0.981 | 0.947 | 0.981 |
| Gis1/Rph1 | B | Rph1 | 3 | 0.921 | 0.973 | 0.940 | 0.980 |
| Gis1/Rph1 | A_B_DBD | Gis1_Rph1_DBD | 4 | 0.932 | 0.950 | 0.953 | 0.966 |
| Gis1/Rph1 | B_A_DBD | Rph1_Gis1_DBD | 4 | 0.905 | 0.870 | 0.934 | 0.914 |
| Pdr8/Yrr1 (foreground) | A | Pdr8 | 2 | 0.736 | 0.924 | 0.736 | 0.924 |
| Pdr8/Yrr1 (foreground) | B | Yrr1 | 2 | 0.729 | 0.957 | 0.729 | 0.957 |
| Pdr8/Yrr1 (foreground) | A_B_DBD | Pdr8_Yrr1_DBD | 3 | 0.681 | 0.838 | 0.742 | 0.873 |
| Pdr8/Yrr1 (foreground) | B_A_DBD | Yrr1_Pdr8_DBD | 4 | 0.599 | 0.814 | 0.672 | 0.871 |
| Vhr1/Vhr2 | A | Vhr1 | 3 | 0.679 | 0.894 | 0.723 | 0.919 |
| Vhr1/Vhr2 | B | Vhr2 | 3 | 0.464 | 0.778 | 0.519 | 0.834 |
| Vhr1/Vhr2 | A_B_DBD | Vhr1_Vhr2_DBD | 3 | 0.937 | 0.970 | 0.951 | 0.978 |
| Vhr1/Vhr2 | B_A_DBD | Vhr2_Vhr1_DBD | 3 | 0.874 | 0.956 | 0.903 | 0.967 |
| Yhp1/Yox1 (foreground) | A | Yhp1 | 2 | 0.860 | 0.975 | 0.860 | 0.975 |
| Yhp1/Yox1 (foreground) | B | Yox1 | 2 | 0.938 | 0.974 | 0.938 | 0.974 |
| Yhp1/Yox1 (foreground) | A_B_DBD | Yhp1_Yox1_DBD | 2 | 0.598 | 0.826 | 0.598 | 0.826 |
| Yhp1/Yox1 (foreground) | B_A_DBD | Yox1_Yhp1_DBD | 3 | 0.889 | 0.922 | 0.914 | 0.941 |

## Correlation of mean tracks across the 4 experiments of each set

Pearson r between replicate-mean tracks, log2(CPM per 100 bp + 1). A = wild-type TF A, B = wild-type TF B, A_B_DBD = A body carrying B's DNA-binding domain, B_A_DBD = B body carrying A's DNA-binding domain.

### Genome-wide (all 100 bp bins, chrM excluded)

| set | A vs B | A vs A_B_DBD | A vs B_A_DBD | B vs A_B_DBD | B vs B_A_DBD | A_B_DBD vs B_A_DBD | n peak bins |
|---|---|---|---|---|---|---|---|
| Dot6/Tod6 (foreground) | 0.871 | 0.945 | 0.789 | 0.887 | 0.827 | 0.805 | 1765 |
| Ecm22/Upc2 | 0.935 | 0.792 | 0.853 | 0.790 | 0.866 | 0.832 | 1856 |
| Fkh1/Fkh2 | 0.864 | 0.962 | 0.830 | 0.874 | 0.958 | 0.838 | 1952 |
| Gis1/Rph1 | 0.749 | 0.962 | 0.806 | 0.723 | 0.948 | 0.812 | 2315 |
| Pdr8/Yrr1 (foreground) | 0.741 | 0.831 | 0.727 | 0.759 | 0.796 | 0.753 | 2616 |
| Vhr1/Vhr2 | 0.542 | 0.879 | 0.674 | 0.617 | 0.783 | 0.775 | 2355 |
| Yhp1/Yox1 (foreground) | 0.911 | 0.795 | 0.901 | 0.785 | 0.954 | 0.797 | 1664 |

### At peaks (union of top-1000 bins of the 4 mean tracks)

| set | A vs B | A vs A_B_DBD | A vs B_A_DBD | B vs A_B_DBD | B vs B_A_DBD | A_B_DBD vs B_A_DBD | n peak bins |
|---|---|---|---|---|---|---|---|
| Dot6/Tod6 (foreground) | 0.248 | 0.935 | 0.285 | 0.251 | 0.803 | 0.272 | 1765 |
| Ecm22/Upc2 | 0.739 | 0.249 | 0.179 | 0.217 | 0.120 | 0.898 | 1856 |
| Fkh1/Fkh2 | -0.015 | 0.914 | -0.046 | 0.067 | 0.949 | -0.013 | 1952 |
| Gis1/Rph1 | -0.349 | 0.934 | -0.174 | -0.466 | 0.884 | -0.181 | 2315 |
| Pdr8/Yrr1 (foreground) | 0.247 | 0.582 | 0.334 | 0.400 | 0.844 | 0.423 | 2616 |
| Vhr1/Vhr2 | -0.244 | 0.896 | -0.117 | -0.100 | 0.918 | 0.028 | 2355 |
| Yhp1/Yox1 (foreground) | 0.614 | 0.650 | 0.621 | 0.439 | 0.905 | 0.543 | 1664 |
