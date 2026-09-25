# Peak-region visualization of GSE179430 ChEC-seq

Environment: `/home/mchikina/work/IDR/.venv` (uv; numpy, pandas, scipy, matplotlib).
Run in order: `01_bin.py` -> `02_select_regions.py` -> `03_extract_windows.py` -> `04_plot.py`.

## Method

1. `01_bin.py`: all 348 per-base coverage files summed into 100 bp bins; also per-bin
   max single-base count and per-sample total counts (`sample_totals.tsv`).
2. `02_select_regions.py`: replicates pooled per experiment (138 experiments), CPM.
   Score per bin per experiment = log2((CPM_300bp-smoothed + 1) / (genome median + 1)).
   Excluded: chrM, rDNA (chrXII:451-469 kb), 5 kb at chromosome ends, experiments with
   < 1M pooled reads, bins where one base holds > 30% of bin counts in the nominating experiment.
   Greedy pick of the 10 top-scoring bins, >= 5 kb apart, each nominated by a different
   experiment. Windows = 2 kb centred on the peak bin. Output: `regions.tsv`.
3. `03_extract_windows.py`: base-resolution counts for all 348 samples in the 10 windows
   (`windows.npz`).
4. `04_plot.py`: figures below and `peak_height_by_experiment.tsv` (max CPM per 25 bp,
   per experiment per region).

## Figures

- `fig_heatmap_10regions.{png,pdf}`: rows = 138 experiments (replicate-pooled; n = replicates),
  grouped by paralog pair (thin separators), ordered wild type / paralog deletion / K. lactis
  ortholog / DBD swap. Columns = 10 regions, 2 kb each at 25 bp resolution, x axis relative
  to the window centre. Colour = CPM per 25 bp, log scale, shared across all panels.
- `fig_lines_10regions.{png,pdf}`: one panel per region. Every experiment as a thin gray
  line; the three experiments with the highest peak in that window coloured and labelled.
  y = CPM per bp, 10 bp running mean; x = genomic coordinate (sacCer3).

## Regions (`regions.tsv`)

| region | locus | nominating experiment |
|---|---|---|
| R1 | chrIII:77,348-79,348 | Stp1_dSTP2 |
| R2 | chrIV:550,128-552,128 | Stp1 |
| R3 | chrVII:270,960-272,960 | Met31_dMET32 |
| R4 | chrXV:618,381-620,381 | Pdr1 |
| R5 | chrVII:609,060-611,060 | Upc2_dECM22 |
| R6 | chrXV:835,781-837,781 | Oaf1_lactis |
| R7 | chrXII:575,522-577,522 | Yap5 |
| R8 | chrVII:986,360-988,360 | Hal9_dTBS1 |
| R9 | chrIII:304,648-306,648 | Tbs1_dHAL9 |
| R10 | chrV:538,295-540,295 | Yap5_dYAP7 |
