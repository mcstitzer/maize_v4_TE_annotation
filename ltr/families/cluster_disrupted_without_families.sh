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

GENOMENAME=B73V4.both_pseudo_AND_unplaced


COMBOLTRFILE=../B73V4.both_pseudo_AND_unplaced.ltrdigest_5ltr.fas
COMBOBASE=$( basename $COMBOLTRFILE .fas )
LTRFILE=../ltrdigest_fastas/B73V4.pseudomolecule.ltrdigest.fastagrab_5ltr.fas
LTRBASE=$( basename $LTRFILE .fas )
UNPLACEDLTR=../ltrdigest_unplaced/B73V4.unplaced.ltrdigest_5ltr.fas
UNPLACEDLTRBASE=$( basename $UNPLACEDLTR .fas )

## make the index for everything first
for i in ../mask_subtract/subtract*/ltrdigest/*_5ltr.fas
do
	samtools faidx $i
done

for i in ../mask_subtract_unplaced/subtract*/ltrdigest/*_5ltr.fas
do
	samtools faidx $i
done

## make sure this is cleared out first so i can append
if [ -f ${GENOMENAME}.nofamily.5ltr.fa ]; then
 rm ${GENOMENAME}.nofamily.5ltr.fa
fi


### don't know how to put the genome name into this quote disaster
### get the chromosome 
awk '{system("samtools faidx ../mask_subtract/subtract"$1"/ltrdigest/B73V4.pseudomolecule.subtract"$1".ltrdigest_5ltr.fas "$2)}' ../copies_without_family_chromosome.txt >> ${GENOMENAME}.nofamily.5ltr.fa
## get the unplaced contigs
awk '{system("samtools faidx ../mask_subtract_unplaced/subtract"$1"/ltrdigest/B73V4.unplaced.subtract"$1".ltrdigest_5ltr.fas "$2)}' ../copies_without_family_unplaced.txt >> ${GENOMENAME}.nofamily.5ltr.fa


${USEARCH} -allpairs_global ${GENOMENAME}.nofamily.5ltr.fa -blast6out ${GENOMENAME}.nofamily.5ltr.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

### so this doesn't get messy, I'll name these families differently as DLTR prefix
${SILIX} ${GENOMENAME}.nofamily.5ltr.fa ${GENOMENAME}.nofamily.5ltr.allvall.8080.out -f DLTR -i 0.8 -r 0.8 --net > ${GENOMENAME}.nofamily.5ltr.8080.fnodes


