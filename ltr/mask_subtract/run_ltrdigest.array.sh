#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/run_ltrdigest_nest-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/run_ltrdigest_nest-stderr-%j.txt
#SBATCH -J ltrdigest
set -e
set -u

export PATH=$PATH:/home/mstitzer/software/bin

### genome tools path
GENOMETOOLS=/home/mstitzer/software/genometools-1.5.1/bin/gt
#GENOMETOOLS=/home/mstitzer/software/genometools-1.5.7/bin/gt
GENOMEBASE=B73V4.pseudomolecule


i=$SLURM_ARRAY_TASK_ID
GENOME=${GENOMEBASE}.subtract${i}

GENOMEFASTA=${GENOME}.fa

MEMLIM=96GB
CPU=16



###################
## run ltrdigest ##
###################

mkdir -p subtract${i}/ltrdigest


echo $GENOME.ltrdigest
echo $GENOME.ltrharvest.sorted.gff3
echo $GENOME
echo $GENOME.ltrdigest.gff3

#$GENOMETOOLS -j $CPU ltrdigest -outfileprefix hardmask${i}/ltrdigest/$GENOME.ltrdigest -v yes -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- hardmask${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME > hardmask${i}/$GENOME.ltrdigest.gff3
## troubleshooting thinking multiple processors weren't working, but looks like it's storing stuff in /tmp on the compute node

### fix the .des files that look funny
#gt encseq encode -des yes -dna yes B73.Mhap2.quiver.subtract13.fa
#$GENOMETOOLS encseq encode -des yes -ssp no -sds no -md5 no -dna yes -indexname $GENOME $GENOMEFASTA

#$GENOMETOOLS suffixerator -db ${GENOMEFASTA} -indexname ${GENOME}.try -des -ssp no -sds no -md5 no -dna -memlimit $MEMLIM

$GENOMETOOLS -j $CPU ltrdigest -outfileprefix subtract${i}/ltrdigest/$GENOME.ltrdigest -trnas ../eukaryotic-tRNAs.fa -hmms ../gydb_hmms/GyDB_collection/profiles/*.hmm -- subtract${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME > subtract${i}/$GENOME.ltrdigest.gff3

## troubleshooting further -- most basic command
##### I DO NOT KNOW WHY BUT I HAVE TO USE genometools-1.5.1 to get the hmm domains to work. I'm giving up on this for now. 

#$GENOMETOOLS -j $CPU ltrdigest -outfileprefix digesttest -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -trnas eukaryotic-tRNAs.fa hardmask${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME
#$GENOMETOOLS -j $CPU ltrdigest -outfileprefix digesttest -hmms all_gydb_profiles.hmmer3b.hmm -trnas eukaryotic-tRNAs.fa hardmask${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME

