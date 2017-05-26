#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract
#SBATCH -o /home/mstitzer/projects/agpv4_te_annotation/slurm-log/subtractgff_nest-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/agpv4_te_annotation/slurm-log/subtractgff_nest-stderr-%j.txt
#SBATCH -J subtractgff
set -e
set -u


export PATH=$PATH:/home/mstitzer/software/bin



i=$SLURM_ARRAY_TASK_ID



#### will generate gffs and gl RDS file for each array
Rscript switch_subtracted_gffs_to_genomic_coordinates.R $i


