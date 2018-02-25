#!/bin/bash -login
#SBATCH -o ltrharvest-stdout-%j.txt
#SBATCH -e ltrharvest-stderr-%j.txt
#SBATCH -J ltrharvest
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=96G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all

set -e
set -u

### genome tools path
module load genometools/1.5.6
module load hmmer/3.1b2
#GENOMETOOLS=/home/mstitzer/software/genometools-1.5.7/bin/gt
GENOMETOOLS=gt

# name the file stem based on suffixator index
GENOME=NIOBT_r1.0
GENOMEFASTA=~/shared/Nobtusifolia/Genome_Files/${GENOME}.fasta

MEMLIM=96GB
CPU=$SLURM_NTASKS

###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################

#$GENOMETOOLS suffixerator -db $GENOMEFASTA -indexname $GENOME -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM

#####################
## Run LTR harvest ##
#####################

#mkdir -p outinner

## all defaults except for maxdistltr (default 15000)
#$GENOMETOOLS ltrharvest -index $GENOME -gff3 $GENOME.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr 20000 -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner outinner/${GENOME}.ltrharvest.outinner.fa -out ${GENOME}.ltrharvest.fa > ${GENOME}.ltrharvest.out

#$GENOMETOOLS gff3 -sort $GENOME.ltrharvest.gff3 > $GENOME.ltrharvest.sorted.gff3

###################
## run ltrdigest ##
###################

#mkdir -p ltrdigest

$GENOMETOOLS -j 16 ltrdigest -outfileprefix ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- $GENOME.ltrharvest.sorted.gff3 $GENOME > $GENOME.ltrdigest.gff3


