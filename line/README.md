## LINE Background
LINEs are a form of nonLTR retrotransposon, which typically contains two ORFs - one for RNA binding, and another for reverse transcriptase and RNaseH. 
Although not every transposition of a LINE generates a TSD, we only annotate copies that have resulted from complete retrotransposition and repair of the target site, not those where TPRT aborts prematurely.


## How LINEs were Identified

The script ```target_mTEA_array.LINE.sh``` takes each LINE reference from the MTEC consensuses, and runs mTEA.
To rerun, need to check paths to all executables, including perl libraries needed for mTEA. 

Each LINE copy with a TSD retained its family assigned by MTEC.

LINEs were concatenated as in ```combine_LINEs.sh```
