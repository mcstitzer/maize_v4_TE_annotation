#!/bin/bash -login
#SBATCH -o run_ltrharvest_nest-stdout-%j.txt
#SBATCH -e run_ltrharvest_nest-stderr-%j.txt
#SBATCH -J ltr_nest_subtract
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=96G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
set -e
set -u

### update March 22, 2016
###   was not removing the copy of the TSD generated via transposition previously. added before and after lines to grep to get the first TSD and the LTR TE itself, which will be removed from the genomic sequence for the next round of search.



MEMLIM=96GB

### load packages
module load genometools/1.5.6
module load bedtools/2.25.0
module load samtools/1.6
GENOMETOOLS=gt

i=1
GENOMEBASE=NIOBT_r1.0

## the first one is different, becuase need to set up subtract directory structure.

python convert_ltrharvest_seq_gff_to_contignames.py ../${GENOMEBASE}.ltrharvest.gff3 > ${GENOMEBASE}.ltrharvest.contignames.gff3
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" ${GENOMEBASE}.ltrharvest.contignames.gff3 | sed -n '1~2p' > ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx ~/shared/Nobtusifolia/Genome_Files/${GENOMEBASE}.fasta
awk -v OFS='\t' '{print $1, $2}' ~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta.fai > ~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta.2.fai
bedtools complement -i ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ~/shared/Nobtusifolia/Genome_Files/${GENOMEBASE}.fasta.2.fai > ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
#fix the gff with 0 start and end lines
grep -Pv "scaffold\d*\t0\t0" NIOBT_r1.0.ltrharvest.contignames.NOTltrretrotransposon.gff3 >NIOBT_r1.0.ltrharvest.contignames.NOTltrretrotransposon.2.gff3
bedtools getfasta -fi ~/shared/Nobtusifolia/Genome_Files/${GENOMEBASE}.fasta -bed ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 -fo ${GENOMEBASE}.subtract1.fa

## concatenate the entries by chromosome
python collapse_chromosomes.py ${GENOMEBASE}.subtract1.fa > ${GENOMEBASE}.temp
mv ${GENOMEBASE}.temp ${GENOMEBASE}.subtract1.fa

### index this fasta
$GENOMETOOLS suffixerator -db ${GENOMEBASE}.subtract1.fa -indexname ${GENOMEBASE}.subtract1 -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
mkdir -p subtract1
$GENOMETOOLS ltrharvest -index ${GENOMEBASE}.subtract1 -gff3 subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr 21000 -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner subtract1/${GENOMEBASE}.subtract1.ltrharvest.outinner.fa -out subtract1/${GENOMEBASE}.subtract1.ltrharvest.fa > subtract1/${GENOMEBASE}.subtract1.ltrharvest.out
$GENOMETOOLS gff3 -sort subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 > subtract1/${GENOMEBASE}.subtract1.ltrharvest.sorted.gff3
#gt gff3 -sort hardmask1/${GENOMEBASE}.hardmask1.ltrharvest.gff3 > hardmask1/${GENOMEBASE}.hardmask1.ltrharvest.sorted.gff3




### go until there are no more LTR TEs in this nested form 
#while [ grep -c ltr_retrotransposon ${GENOMEBASE}.hardmask${i}.ltrharvest.gff3 -gt 0 ]
while [ $i -le 5 ] #was 100
do

# name the file stem based on suffixator index

OLDINDEX=$i
GENOME=${GENOMEBASE}.subtract${i}
i=$(( $i + 1 ))
NEWINDEX=$i
NEWGENOME=${GENOMEBASE}.subtract${i}

GENOMEFASTA=${GENOME}.fa
NEWGENOMEFASTA=${NEWGENOME}.fa

MEMLIM=96GB
CPU=$SLURM_NTASKS

###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################


## switch genome tools back to their real contig names
python convert_ltrharvest_seq_gff_to_contignames.py subtract${OLDINDEX}/${GENOME}.ltrharvest.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3
## only get the LTR_retrotransposon records to keep Rscript from repeating computation
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx ${GENOMEFASTA}
awk -v OFS='\t' '{print $1, $2}' ${GENOMEFASTA}.fai > ${GENOMEFASTA}.1.fai
mv ${GENOMEFASTA}.1.fai ${GENOMEFASTA}.fai 
### run the rscript: read in gff, read in RDS of genomelist, update genomelist with updatePos, write gff with changed positions


## find regions not covered by TEs
bedtools complement -i subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ${GENOME}.fa.fai > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
grep -Pv "scaffold\d*\t0\t0" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3
mv subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3 subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3
bedtools getfasta -fi $GENOMEFASTA -bed subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo $NEWGENOMEFASTA

## concatenate the entries by chromosome
python collapse_chromosomes.py $NEWGENOMEFASTA > ${NEWGENOMEFASTA}.temp
mv ${NEWGENOMEFASTA}.temp $NEWGENOMEFASTA

### index this fasta
$GENOMETOOLS suffixerator -db ${NEWGENOMEFASTA} -indexname ${NEWGENOME} -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM


#####################
## Run LTR harvest ##
#####################

mkdir -p subtract${NEWINDEX}
### allow extra 1kb for each iteration, because we miss insertions that are not structural
MAXLEN=$(($i * 1000 + 20000))    ### so 20kb for the first hardmask, plus the additional 1kb per round
## all defaults except for maxdistltr (default 15000)
$GENOMETOOLS ltrharvest -index ${NEWGENOME} -gff3 subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr $MAXLEN -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.outinner.fa -out subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.fa > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.out

$GENOMETOOLS gff3 -sort subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.sorted.gff3


done


## can then run all the ltrdigest in array form on all the gffs
###################
## run ltrdigest ##
###################

#mkdir -p hardmask5/ltrdigest

#$GENOMETOOLS -j 16 ltrdigest -outfileprefix hardmask5/ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- hardmask5/$GENOME.ltrharvest.sorted.gff3 $GENOME > hardmask5/$GENOME.ltrdigest.gff3

