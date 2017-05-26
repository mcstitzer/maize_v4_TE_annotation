#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/run_ltrharvest_nest-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/run_ltrharvest_nest-stderr-%j.txt
#SBATCH -J ltr_nest_subtract
set -e
set -u



### update March 22, 2016
###   was not removing the copy of the TSD generated via transposition previously. added before and after lines to grep to get the first TSD and the LTR TE itself, which will be removed from the genomic sequence for the next round of search.

export PATH=$PATH:/home/mstitzer/software/bin

MEMLIM=96GB

### genome tools path
GENOMETOOLS=/home/mstitzer/software/genometools-1.5.7/bin/gt

i=1
GENOMEBASE=B73V4.pseudomolecule


## the first one is different, becuase need to set up subtract directory structure.

python convert_ltrharvest_seq_gff_to_contignames.tabsep.py ../${GENOMEBASE}.ltrdigest.gff3 > ${GENOMEBASE}.ltrdigest.contignames.gff3


while [ $i -le 80 ]
do


OLDINDEX=$i
GENOME=${GENOMEBASE}.subtract${i}
i=$(( $i + 1 ))
NEWINDEX=$i
NEWGENOME=${GENOMEBASE}.subtract${i}

GENOMEFASTA=${GENOME}.fa
NEWGENOMEFASTA=${NEWGENOME}.fa

MEMLIM=96GB
CPU=16

python convert_ltrharvest_seq_gff_to_contignames.tabsep.py subtract${OLDINDEX}/${GENOME}.ltrdigest.gff3 > subtract${OLDINDEX}/${GENOME}.ltrdigest.contignames.gff3


done


