#!/bin/bash

## get tRNA database - all eukaryotes
wget http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
gzip -d eukaryotic-tRNAs.fa.gz
#rm eukaryotic-tRNAs.fa.gz

## get gypdb HMMs
wget http://gydb.org/gydbModules/collection/collection/db/GyDB_collection.zip
mkdir gydb_hmms
unzip GyDB_collection.zip -d gydb_hmms
rm GyDB_collection.zip

### there's a problem with one of them - ty1/copia is the name, and the / is interpreted weird when put into paths
## AP_ty1copia.hmm
### here's the fix, looks funny because wanting to replace a /:
sed -i "s#ty1/copia#ty1-copia#g" gydb_hmms/GyDB_collection/profiles/AP_ty1copia.hmm 

### also, galadriel doesn't have INT_ prefaced in the hmm name. This complicates using the protein domain.

sed -i "s/galadriel/INT_galadriel/" gydb_hmms/GyDB_collection/profiles/INT_galadriel.hmm

