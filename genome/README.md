# S. cerevisiae S288C reference genome, R64 (sacCer3)

Downloaded from UCSC on 2026-09-25:
https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/

UCSC sacCer3 = SGD R64-1-1 = NCBI GCA_000146045.2. This is the assembly the
GSE179430 ChEC-seq coverage files were aligned to (GEO: "S2988C_R64",
"GCA_000146045.2 (R64)"). Genome sequence is unchanged in later R64 annotation
releases (R64-2-1 onward).

Files:
- sacCer3.fa                    uncompressed FASTA, 17 sequences (chrI..chrXVI, chrM), soft-masked (lowercase repeats)
- sacCer3.fa.fai                samtools faidx index
- sacCer3.chrom.sizes           UCSC chromosome sizes (sorted by size, as distributed)
- sacCer3.chrom.sizes.geo_order same, reordered chrI..chrXVI, chrM: the concatenation order of the GSE179430 *.out.txt.gz files

Total length 12,157,105 bp, equal to the number of lines in each GSE179430 coverage file.
