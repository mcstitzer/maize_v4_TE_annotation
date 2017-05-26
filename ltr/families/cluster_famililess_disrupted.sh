#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/families
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/run_ltrharvest-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/run_ltrharvest-stderr-%j.txt
#SBATCH -J subtracted_fam
set -e
set -u


CPU=4


GENOMENAME=B73V4.unplaced
#USEARCH=~/software/usearch8.0.1623_i86linux32 
USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
LTRFILE=../B73V4.both_pseudo_AND_unplaced.ltrdigest_5ltr.fas

i=$SLURM_ARRAY_TASK_ID

#$USEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db ../ltrdigest/B73.Mhap2.quiver.ltrdigest_complete.fas -id 0.9 -blast6out B73.Mhap2.quiver.ltrdigest.TEDB.searchglobal.out -strand both -query_cov 0.9 -target_cov 0.9 --threads 24

mkdir -p subtracted
$USEARCH --usearch_global ../mask_subtract_unplaced/subtract${i}/ltrdigest/${GENOMENAME}.subtract${i}.ltrdigest_5ltr.fas -db $LTRFILE -id 0.8 -blast6out subtracted/${GENOMENAME}.ltrdigest.subtract${i}.ORIG5pLTR.id80.cov80.out -strand both -query_cov 0.8 -target_cov 0.8 -top_hits_only --threads $CPU
#$USEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db ../ltrdigest/B73.Mhap2.quiver.ltrdigest_complete.fas -id 0.9 -blast6out B73.Mhap2.quiver.ltrdigest.TEDB.searchglobal.query90.out -strand both -query_cov 0.9 --threads 16
#$USEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db $LTRFILE -id 0.9 -blast6out B73.Mhap2.quiver.ltrdigest.TEDB.5ltr.out -strand both -target_cov 0.9 --threads 16

#makeblastdb -in ../ltrdigest/B73.Mhap2.quiver.ltrdigest_complete.fas -dbtype nucl
#blastn -query TE_12-Feb-2015_15-35.fa -db ../ltrdigest/B73.Mhap2.quiver.ltrdigest_complete.fas -outfmt 6 -max_target_seqs 1 -perc_identity 95 -word_size 100 -out B73.Mhap2.quiver.ltrdigest.TEDB.blast.out -num_threads 16

