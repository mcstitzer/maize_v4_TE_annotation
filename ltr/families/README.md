## clustering LTR TEs into families

**Note** this is not what I would do again, but scripts here to replicate what I did.

To recreate again, I would concatenate all 5' LTRs of all TE copies, and run the clustering on that file. 

- ```cluster_families_LTR.sh``` is a script to run an all-v-all homology search and cluster results.

- ```cluster_mtec.sh``` identifies the MTEC consensus that best matches each family.

Download the [Maize TE Consortium Consensuses](http://maizetedb.org) [here](http://maizetedb.org/~maize/TE_12-Feb-2015_15-35.fa)
