#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/families
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/cluster_families-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/cluster_families-stderr-%j.txt
#SBATCH -J cluster_LTR
set -e
set -u


CPU=24

#USEARCH=~/software/usearch8.0.1623_i86linux32 
USEARCH=/home/mstitzer/software/vsearch/bin/vsearch
SILIX=~/software/bin/silix
COMBOLTRFILE=../B73V4.both_pseudo_AND_unplaced.ltrdigest_5ltr.fas
COMBOBASE=$( basename $COMBOLTRFILE .fas )
LTRFILE=../ltrdigest_fastas/B73V4.pseudomolecule.ltrdigest.fastagrab_5ltr.fas
LTRBASE=$( basename $LTRFILE .fas )
UNPLACEDLTR=../ltrdigest_unplaced/B73V4.unplaced.ltrdigest_5ltr.fas
UNPLACEDLTRBASE=$( basename $UNPLACEDLTR .fas )



#${USEARCH} -allpairs_global ${LTRFILE} -blast6out ${LTRBASE}.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

#${USEARCH} -allpairs_global ${UNPLACEDLTR} -blast6out ${UNPLACEDLTRBASE}.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

#${USEARCH} --usearch_global $LTRFILE -db $UNPLACEDLTR -blast6out ${LTRBASE}.v.${UNPLACEDLTRBASE}.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 -strand both --threads $CPU

#${USEARCH} --usearch_global $UNPLACEDLTR -db $LTRFILE -blast6out ${UNPLACEDLTRBASE}.v.${LTRBASE}.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 -strand both --threads $CPU

#cat ${LTRBASE}.allvall.8080.out ${UNPLACEDLTRBASE}.allvall.8080.out ${LTRBASE}.v.${UNPLACEDLTRBASE}.8080.out ${UNPLACEDLTRBASE}.v.${LTRBASE}.8080.out > ${COMBOBASE}.combined.allvall.8080.out


${SILIX} ${COMBOLTRFILE} ${COMBOBASE}.combined.allvall.8080.out -f LTR -i 0.8 -r 0.8 --net > ${COMBOBASE}.combined.8080.fnodes


