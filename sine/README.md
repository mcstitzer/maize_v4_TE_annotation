## SINE Background
SINEs are transcribed by RNA polymerase III, and are derived from one of three classes of Pol III–transcribed molecules (tRNA, 7SL, 5s rRNA).
While animal SINEs from all three classes are known, plant SINEs are exclusively derived from tRNA. 
To find SINEs, I apply the implementation of [Wenke et al. 2011](http://www.plantcell.org/content/early/2011/09/07/tpc.111.088682) in [SINE-finder](http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt), which searches for tRNA-derived SINEs containing RNA polymerase III A and B boxes near the polyA tail. 
The defaults are that A and B box consensus nucleotide sequences are RVTGG and GTTCRA, there is a 25–50 bp spacer between the A and B boxes, and there is a spacer of 20–500 bp between the B box and polyA tail.

## How SINEs were Identified

The script ```run_sines.sh``` will download sine_finder, run sine_finder, parse results, assign to families, and output a GFF.
To rerun, need to check paths at the beginning of the file, where things like the genome name and paths to executables for SILIX and VSEARCH are located.

Each candidate SINE was is clustered using VSEARCH and silix, to characterize families.
These families are clustered with [Maize TE Consortium](http://www.maizetedb.org) exemplars to faciliatate comparison between genome versions.

The script ```slurm_sinefinder.sh``` is what I used to submit on our cluster.
