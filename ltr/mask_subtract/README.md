## remove existing LTR TEs and find sandwiched LTR TEs

- the script ```run_ltrharvest.nestedloop.sh``` removes previous TEs, makes a subtracted genome, and reruns, making use of ```collapse_chromosomes.py``` and ```convert_ltrharvest_seq_gff_to_contignames.py```.

- On our cluster running SLURM, ```run_ltrdigest.array.sh``` submits an array job to run ltrdigest on the TE models from each subtracted round.

- Then, we need to convert these subtracted coordinates back to real genomic coordinates and chromosome names, using ```submit_subtraction_array.sh```. This makes use of ```switch_subtracted_gffs_to_genomic_coordinates.R```

