#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/maize_v4_TE_annotation/sine
#SBATCH -o /home/mstitzer/projects/maize_v4_TE_annotation/slurm-log/sine_finder-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/maize_v4_TE_annotation/slurm-log/sine_finder-stderr-%j.txt
#SBATCH -J sine_finder
set -e
set -u



chmod 755 run_sines.sh 
./run_sines.sh
