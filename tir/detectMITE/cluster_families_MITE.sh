#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/detectMITE
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/cluster_families-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/cluster_families-stderr-%j.txt
#SBATCH -J cluster_MITE
set -e
set -u


CPU=24

#USEARCH=~/software/usearch8.0.1623_i86linux32 
USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
MITEFILE=B73V4.both_pseudo_AND_unplaced.miteSet.fasta
MITEBASE=$( basename $MITEFILE .fasta )

MITEOUT=${MITEBASE}.nodash.fasta

if [ ! -f $MITEOUT ];
then
        grep -v -- "----" ${MITEFILE} > ${MITEOUT}
fi


#${USEARCH} -allpairs_global ${MITEOUT} -blast6out ${MITEBASE}.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

${SILIX} ${MITEOUT} ${MITEBASE}.allvall.8080.out -f MITE -i 0.8 -r 0.8 --net > ${MITEBASE}.8080.fnodes


