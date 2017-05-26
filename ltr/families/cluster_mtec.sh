#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/families
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/mtec_fam-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/mtec_fam-stderr-%j.txt
#SBATCH -J mtec_fam
set -e
set -u

CPU=24

GENOMENAME=B73V4.both_pseudo_AND_unplaced

USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
TEFILE=../${GENOMENAME}.ltrdigest_complete.fas


$USEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db $TEFILE -id 0.8 -query_cov 0.8 -target_cov 0.8 -blast6out ${GENOMENAME}.TEDB.8080.searchglobal.toponly.out -strand both -top_hits_only --threads $CPU

