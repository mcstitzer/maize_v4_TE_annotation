#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/tir/
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/cluster_families-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/cluster_families-stderr-%j.txt
#SBATCH -J cluster_tir
set -e
set -u


CPU=24

#USEARCH=~/software/usearch8.0.1623_i86linux32 
USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
TIRFILE=B73V4.both_pseudo_AND_unplaced.tirmite.fa
TIRBASE=$( basename $TIRFILE .fa )

${USEARCH} -allpairs_global ${TIRFILE} -blast6out ${TIRBASE}.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

${SILIX} ${TIRFILE} ${TIRBASE}.allvall.8080.out -f TIR -i 0.8 -r 0.8 --net > ${TIRBASE}.8080.fnodes


