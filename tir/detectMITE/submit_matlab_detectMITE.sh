#!/bin/bash -login
#SBATCH -D /home/mstitzer/software/detectMITE
#SBATCH -o /home/mstitzer/software/detectMITE/slurm-log/detectMITE-stdout-%j.txt
#SBATCH -e /home/mstitzer/software/detectMITE/slurm-log/detectMITE-stderr-%j.txt
#SBATCH -J detectMITE
set -e
set -u

### 12 cpu, 36000 memory
### this expects https://sourceforge.net/projects/detectmite/files/detectMITE.20170425.tar.gz
### with cd-hit executable inside (detectMITE/cd-hit/cd-hit-est)


module load matlab
#module load matlab matlab/7.11
export LD_LIBRARY_PATH=~/software/detectMITE/bin/glnxa64/:~/software/detectMITE/runtime/

#matlab -nodisplay -nosplash -r "tic;do_MITE_detectionRestart('B73.Mhap2.quiver.fasta','-genome','Mhap2restart','-cpu',12);runtime=toc;quit"
#matlab -nodisplay -nosplash -r "tic;path(path, 'bioinfo/+bioinfo/+bioinfoprivate', 'bioinfo/', 'bioinfo/+bioinfo/');do_MITE_detectionRestart('B73.Mhap2.quiver.fasta','-genome','Mhap2restart','-cpu',12);runtime=toc;quit"
#matlab -nodisplay -nosplash -r "tic;path(path, 'bioinfo/');do_MITE_detectionRestart('B73.Mhap2.quiver.fasta','-genome','Mhap2restart','-cpu',12);runtime=toc;quit"


#matlab -nodisplay -nosplash -r "tic;do_MITE_detectionRestart('utg10007.fa','-genome','utg10007','-cpu',12);runtime=toc;quit"

#matlab -nodisplay -nosplash -r "tic;do_MITE_detectionRestart('B73.Mhap2.quiver.fasta','-genome','Mhap2restart','-cpu',12);runtime=toc;quit"

matlab -nodisplay -nosplash -r "tic;do_MITE_detection('B73V4.both_pseudo_AND_unplaced.fa','-genome','B73V4.both_pseudo_AND_unplaced','-cpu',12);runtime=toc;quit"

#matlab -nodisplay -nosplash -r "tic;do_MITE_detectionMCS('B73.Mhap2.quiver.fasta','-genome','Mhap2_10kb','-cpu',12,'-mite_maximum_length',10000);runtime=toc;quit"

############  GENOME
#matlab -nodisplay -nosplash -r "tic;do_MITE_detectionMCS('B73V4.both_pseudo_AND_unplaced.fa','-genome','B73V4.both_pseudo_AND_unplaced','-cpu',12);runtime=toc;quit"

#matlab -nodisplay -nosplash -r "tic;do_MITE_detectionMCS('B73V4.both_pseudo_AND_unplaced.fa','-genome','B73V4.both_pseudo_AND_unplaced_10kb','-cpu',12,'-mite_maximum_length',10000);runtime=toc;quit"
