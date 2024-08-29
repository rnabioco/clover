## Files

yeast/grande/grande.bam Grandeyeast004_20240105_1222_P2S-00519-B_PAS98845_231d7d7e.rna004_130bps_sup@v5.0.0.bwa.bam
yeast/petite/petite.bam Petityeast004_20240105_1222_P2S-00519-A_PAQ47538_49891fce.rna004_130bps_sup@v5.0.0.bwa.bam
yeast/trna-ref.fa.gz    sacCer-mito-and-nuclear-tRNAs.fa.gz
yeast/trna-mods.tsv.gz  scerevisiae.mods.seq2structure.tsv.gz

## bcerr

detectrms --bam petite.bam --fasta trna-ref.fa -o petite.bcerr.tsv
detectrms --bam grande.bam --fasta trna-ref.fa -o grande.bcerr.tsv
