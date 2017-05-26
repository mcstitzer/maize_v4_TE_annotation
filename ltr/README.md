
## To predict structural LTRs:

- download tRNA and GyDb HMMs using ```get_tRNA_hmm_dbs.sh```, which are needed for ltrdigest

- ```run_ltrharvest.sh``` runs [ltrharvest](http://www.zbh.uni-hamburg.de/fileadmin/gi/LTRharvest/ltrharvestman.pdf) and [ltrdigest](http://www.zbh.uni-hamburg.de/fileadmin/gi/LTRdigest/ltrdigestman.pdf) on the genome

- but LTR TEs are nested, so we need to remove these copies and rerun. This is done in ```mask_subtract```


## To cluster into families

We use the 808080 rule ([Wicker et al., 2007](http://www.nature.com/nrg/journal/v8/n12/full/nrg2165.html)) to cluster. 

In the ```families``` directory:

- `cluster_families_LTR.sh` runs the initial allvsall blast of 5' LTRs and then goes through silix clustering of sequences


sed -i 's/=//g' B73.Mhap2.quiver.ltrdigest_5ltr.fas


	
