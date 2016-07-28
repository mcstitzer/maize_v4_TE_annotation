#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/maize_v4_TE_annotation/helitron
#SBATCH -o /home/mstitzer/projects/maize_v4_TE_annotation/slurm-log/helitron-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/maize_v4_TE_annotation/slurm-log/helitron-stderr-%j.txt
#SBATCH -J helitron
set -e
set -u

chmod 755 run_helitron_scanner.sh

srun ./run_helitron_scanner.sh

