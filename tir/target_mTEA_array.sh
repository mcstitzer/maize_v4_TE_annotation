#!/bin/bash
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/tir/target_mTEA
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/target_dna-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/target_dna-stderr-%j.txt
#SBATCH -J target_dna
#SBATCH -p serial
#SBATCH --ntasks-per-node 16
#SBATCH --mem 24000

module load blast/2.2.26
export PATH=$PATH:/home/mstitzer/software/bin
export PATH=$PATH:/home/mstitzer/projects/agpv4_te_annotation/tir/mTEA/lib/blogo/
#export $PERL5LIB
export PERL5LIB=$PERL5LIB:/home/mstitzer/projects/agpv4_te_annotation/tir/mTEA/lib/blogo/


FILES=($(ls -1 ~/te_reference_fasta/individual_fasta/DTC_ZM*))

### set up a file name for each individual TE fasta
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}
echo $FILENAME

### RUN TARGeT
FILE=$(basename "$FILENAME")
echo $FILE
DNATE="${FILE%.*}"
echo $DNATE

if [ ! -f ${DNATE}.tir.fa.tab ]
then

mkdir -p $DNATE
python ~/software/TARGeT/target.py -q $FILENAME -t nucl -o $DNATE -i s -P 16 -S PHI -DB -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.3 ../../B73.Mhap2.quiver.fasta ${DNATE}_target


### Convert names of flanking fasta file
python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank_adj > ${DNATE}.flank.fa

### RUN mTEA
## CACTA : -c NNN -t CACTNNNNNNNNN -d 3 -s 0
perl ../mTEA/other-scripts/id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNN -t CACTNNNNNNNNN -d 3 -s 0


## hAT : -c NNNNNNNN -t NNNNNNNNNNN -d 2 -s 1
#perl mTEA/other-scripts/id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNNNNNNN -t NNNNNNNNNNN -d 2 -s 1

## Mutator : -c NNNNNNNNNNN -t NNNNNNNNNNN -d 4 -s 1
#perl mTEA/other-scripts/id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNN -t CACTNNNNNNNNN -d 3 -s 0

## Tc1/Mariner : -c TA -t NNNNNNNNNNNN -d 2 -s 0
#perl mTEA/other-scripts/id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c TA -t NNNNNNNNNNNN -d 2 -s 0

## PIF/Harbinger : -c TNN -t NNNNNNNNNNNNNN -d 3 -s 0
#perl mTEA/other-scripts/id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c TNN -t NNNNNNNNNNNNNN -d 3 -s 0

fi

### so to submit, grab files and export them
####### DO ALL OF THIS IN THE SHELL YOU'RE ABOUT TO SUBMIT IN
### export FILES=($(ls -1 DTA*))
### NUMFILES=${#FILES[@]}
### ARRAYNUM=$(($NUMFILES - 1)) ## because slurm array is zero based

## if [ $ARRAYNUM -ge 0 ]; then
## sbatch --array=0-$ARRAYNUM tir_search_genome.sh
## fi
