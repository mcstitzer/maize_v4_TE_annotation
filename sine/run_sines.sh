#!/bin/bash -login

GENOME=B73V4
GENOMEFASTA=../B73V4.fa
CPU=16

## path to SILIX for clustering
SILIX=/home/mstitzer/software/bin/silix
## path to VSEARCH for matching sequences
VSEARCH=/home/mstitzer/software/vsearch/bin/vsearch

#### get the SINE-Finder program from Wenke et al. 2011
wget http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt 

#### change name
mv Supplemental_Data_Set_1-sine_finder.txt sine_finder.py

#### run sinefinder
#### I haven't been able to get sine_finder to work with reverse sequences, as it seems to report TSDs wrong on the reverse strand.
####   so I'm only reporting on the forward strand.
### -f both : outputs csv and fasta
python sine_finder.py -T chunkwise -V1 -f both -o F ../${GENOMEFASTA}

#### sine_finder outputs the fasta with the TSD included. I remove these here, so they aren't considered when clustering into families
mv ../${GENOME}-matches.fasta .
mv ../${GENOME}-matches.csv .
python remove_tsd_sinefinder.py ${GENOME}-matches.fasta ${GENOME}-matches.noTSD.fa

#### vsearch to identify homology, silix to cluster
$VSEARCH -allpairs_global ${GENOME}-matches.noTSD.fa -blast6out ${GENOME}-matches.noTSD.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU

# single linkage cluster those that are 80% identical to each other.
$SILIX ${GENOME}-matches.noTSD.fa -f SINE -i 0.8 -r 0.8 > ${GENOME}-matches.noTSD.8080.fnodes

### cluster my families into MTEC TE families
wget http://maizetedb.org/~maize/TE_12-Feb-2015_15-35.fa
$VSEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db ${GENOME}-matches.noTSD.fa -id 0.8 -query_cov 0.8 -target_cov 0.8 -blast6out ${GENOME}-matches.noTSD.TEDB8080.out -strand both -top_hits_only --threads $CPU

### cluster into families and output final gff with this R script
Rscript generate_gff_SINE.R $GENOME


