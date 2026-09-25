# Prediction targets: mean ChEC-seq tracks for complete DBD-swap sets

Script: `05_targets_and_correlations.py` (env `/home/mchikina/work/IDR/.venv`; needs `../peaks/bins100.npz` from `../peaks/01_bin.py`).

## Sets (7 pairs x 4 constructs = 28 experiments)

Pairs with all four of A, B, A body + B DBD (`A_B_DBD`), B body + A DBD (`B_A_DBD`):
Dot6/Tod6, Ecm22/Upc2, Fkh1/Fkh2, Gis1/Rph1, Pdr8/Yrr1, Vhr1/Vhr2, Yhp1/Yox1.
Ace2/Swi5 and Cup9/Tos8 have only one swap direction in GEO and are excluded.
`manifest.tsv`: set, role, expname, replicate count, source files, per-replicate total counts, track path.

## Mean tracks (`mean_tracks/<expname>.mean_cpm.npy`)

- float32, 12,157,105 values = per-base CPM, genome in GEO order (chrI..chrXVI, chrM; see `../../genome/sacCer3.chrom.sizes.geo_order`).
- Each replicate normalised to its own depth (counts / total x 1e6), then averaged with equal weight across replicates.
- 48.6 MB each, 1.3 GB total. Load with `np.load(path, mmap_mode="r")`.

## Correlations

All Pearson r on log2(CPM per 100 bp bin + 1); chrM excluded.
Peaks per set = union over the set's 4 mean tracks of each track's 1000 highest bins
(excluding chrM, rDNA chrXII:451-469 kb, 5 kb chromosome ends, single-base pile-up bins).
Union sizes are 1,664-2,616 bins per set (`n_peak_bins` in `mean_track_correlation.tsv`).

- `replicate_correlation.tsv`: within each experiment, every replicate pair (`rep_i`, `rep_j`), and every
  replicate vs the mean of the remaining replicates (`rep_j = mean_of_others`), genome-wide (`r_genome`)
  and over the set's peak bins (`r_peaks`).
- `mean_track_correlation.tsv`: within each set, all 6 pairs of mean tracks, `r_genome` and `r_peaks`.

## Foreground sets (decision 2026-09-25)

Foreground = sets in which all 6 mean-track cross-peak correlations (`r_peaks`) exceed 0.2:
**Dot6/Tod6, Pdr8/Yrr1, Yhp1/Yox1** (12 experiments, 31 replicates).
Flagged in `manifest.tsv` column `all_peak_r_above_0.2`.

| set | min r_peaks | weakest comparison |
|---|---|---|
| Dot6/Tod6 | 0.25 | A vs B |
| Pdr8/Yrr1 | 0.25 | A vs B |
| Yhp1/Yox1 | 0.44 | B vs A_B_DBD |

Background (kept on disk, not foreground): Ecm22/Upc2 (min 0.12), Fkh1/Fkh2 (-0.05), Vhr1/Vhr2 (-0.24), Gis1/Rph1 (-0.47).

Foreground experiments and mean tracks:

| set | A | B | A_B_DBD | B_A_DBD |
|---|---|---|---|---|
| Dot6/Tod6 | Dot6 (n=4) | Tod6 (n=2) | Dot6_Tod6_DBD (n=2) | Tod6_Dot6_DBD (n=3) |
| Pdr8/Yrr1 | Pdr8 (n=2) | Yrr1 (n=2) | Pdr8_Yrr1_DBD (n=3) | Yrr1_Pdr8_DBD (n=4) |
| Yhp1/Yox1 | Yhp1 (n=2) | Yox1 (n=2) | Yhp1_Yox1_DBD (n=2) | Yox1_Yhp1_DBD (n=3) |

Tracks: `mean_tracks/<expname>.mean_cpm.npy`.
