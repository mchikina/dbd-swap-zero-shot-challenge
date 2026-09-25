# GSE179430 — WGD paralog evolution (ChEC-seq, S. cerevisiae)

Downloaded from NCBI GEO on 2026-09-25.

## Source

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179430
- BioProject: PRJNA743682
- SRA: SRP326882 (raw reads; NOT downloaded here)
- Publication: Gera T, Jonas F, More R, Barkai N. *Evolution of binding preferences
  among whole-genome duplicated transcription factors.* eLife 2022; 11:e73225.
  PMID 35404235, doi:10.7554/eLife.73225
- Contact: Naama Barkai lab, Weizmann Institute
- Public since 2022-01-21; last GEO update 2022-04-20

## Summary (from GEO)

Studying the evolution of binding preferences for 30 paralog TF pairs in
S. cerevisiae. Genome-wide binding of wild-type TFs and orthologs in different
backgrounds was profiled using ChEC-seq.

## Directory layout

```
GSE179430/
  README.md
  GSE179430_RAW.tar                 original GEO archive (469,524,480 bytes, md5 917b9f0d952802b5497d748fb5ccaa5a)
  raw/                              348 files extracted from the tar, one per GSM sample
    GSM5417602_Cup9_rpt1.out.txt.gz
    ...
    GSM5417949_Fkh2_lactis_rpt2.out.txt.gz
  metadata/
    samples.tsv                     one row per sample: gsm, title, expname, strainid, genotype, platform, instrument, file
    filelist.txt                    GEO listing of tar contents with sizes
    GSE179430_family.soft.txt       full GEO SOFT record (series, platforms, all 348 samples)
    GSE179430-GPL19756_series_matrix.txt.gz   GEO series matrix, NextSeq 500 samples (19)
    GSE179430-GPL27812_series_matrix.txt.gz   GEO series matrix, NovaSeq 6000 samples (329)
```

All 348 extracted files match the sizes in `metadata/filelist.txt`.

## Processed data format (`raw/*.out.txt.gz`)

Per GEO: "*.out file contains genome coverage for all 12,157,105 bases in the
yeast genome in one column (chromosomal order 1-16, Mitochondrion)".

- Single column of integers, no header, one line per genomic base.
- Value = number of read 5' ends at that base (ChEC MNase cut sites).
- 12,157,105 lines per file (verified on GSM5417602).
- Genome build: GCA_000146045.2 (R64 / sacCer3). Reference FASTA and chrom.sizes
  are in `../genome/` (see its README). Concatenation order and lengths
  (sum = 12,157,105):

| chr   | length    | first line | last line  |
|-------|-----------|------------|------------|
| I     |   230,218 |          1 |    230,218 |
| II    |   813,184 |    230,219 |  1,043,402 |
| III   |   316,620 |  1,043,403 |  1,360,022 |
| IV    | 1,531,933 |  1,360,023 |  2,891,955 |
| V     |   576,874 |  2,891,956 |  3,468,829 |
| VI    |   270,161 |  3,468,830 |  3,738,990 |
| VII   | 1,090,940 |  3,738,991 |  4,829,930 |
| VIII  |   562,643 |  4,829,931 |  5,392,573 |
| IX    |   439,888 |  5,392,574 |  5,832,461 |
| X     |   745,751 |  5,832,462 |  6,578,212 |
| XI    |   666,816 |  6,578,213 |  7,245,028 |
| XII   | 1,078,177 |  7,245,029 |  8,323,205 |
| XIII  |   924,431 |  8,323,206 |  9,247,636 |
| XIV   |   784,333 |  9,247,637 | 10,031,969 |
| XV    | 1,091,291 | 10,031,970 | 11,123,260 |
| XVI   |   948,066 | 11,123,261 | 12,071,326 |
| M     |    85,779 | 12,071,327 | 12,157,105 |

## Processing pipeline (from GEO sample records)

- Library strategy: ChEC-seq (Zentner et al. 2015); one-tube library construction
  based on Skene et al. 2017 as described in Gera 2021
- Growth: overnight cultures diluted >1:10000, grown to OD 4 in SC + glucose
- Demultiplexing: bcl2fastq
- Adapter-dimer filtering: cutadapt
- Alignment: bowtie2, paired-end, to S288C R64
  (`-p8 --local --very-sensitive --trim-to 30 --dovetail`)
- Coverage: bedtools genomecov (`-5 -fs 1 -d`)
- Instruments: Illumina NovaSeq 6000 (GPL27812, 329 samples), NextSeq 500 (GPL19756, 19 samples)

## Sample naming

Sample titles / file names follow `<expname>_rpt<N>`, where `<expname>` encodes
the construct (see `metadata/samples.tsv` for the exact genotype of each strain).
All strains are BY4741 with the TF fused to MNase at its endogenous locus.

| pattern                 | n   | meaning (from genotype field)                                    | example genotype                          |
|-------------------------|-----|------------------------------------------------------------------|-------------------------------------------|
| `TF`                    | 118 | wild-type TF-MNase                                               | `BY4741 CUP9::Cup9-Mnase`                 |
| `TF_dPARALOG`           | 131 | TF-MNase with its paralog deleted (kanMX)                        | `BY4741 GIS1::Gis1-Mnase RPH1::kanMX`     |
| `TF_lactis`             |  48 | K. lactis ortholog-MNase at the TF locus, S. cerevisiae paralog deleted | `BY4741 OAF1::K.lactis A03443-Mnase PIP2::kanMX` |
| `TF_PARALOG_DBD`        |  51 | TF-MNase with its DNA-binding domain swapped for the paralog's   | `BY4741 YRR1::Yrr1-Pdr8DBD-Mnase`         |

138 distinct experiment names (expname); 1–5 replicates each.

Paralog pairs, as encoded in the `_d<PARALOG>` / `_<PARALOG>_DBD` sample names
(29 pairs): Cup9/Tos8, Skn7/Hms2, Dot6/Tod6, Sok2/Phd1, Smp1/Rlm1, Gzf3/Dal80,
Gis1/Rph1, Met31/Met32, Stp3/Stp4, Stp1/Stp2, Pdr1/Pdr3, Vhr1/Vhr2, Pip2/Oaf1,
Ixr1/Abf2, Ecm22/Upc2, Yrr1/Pdr8, Hal9/Tbs1, Tda9/Rsf2, Yap6/Yap4, Yap5/Yap7,
Yap1/Yap2, Yhp1/Yox1, Sut1/Sut2, Spt23/Mga2, Msn4/Msn2, Mig2/Mig3, Nrg1/Nrg2,
Fkh1/Fkh2, Ace2/Swi5. Usv1 and Rgm1 appear only as wild-type samples; their
pairing is not encoded in the sample names.

Only one member of some pairs has a wild-type sample (e.g. Abf2, Msn2, Mig3,
Nrg1, Nrg2, Yap1, Yap2 have no plain `TF_rptN` sample).

## Not downloaded

- Raw FASTQ / SRA runs (SRP326882).
