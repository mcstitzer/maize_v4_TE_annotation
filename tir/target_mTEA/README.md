## Terminal Inverted Repeat (TIR) Background
There are 5 superfamilies of TIR elements in plants: hAT (DTA), CACTA (DTC), pIF/Harbinger (DTH), Mutator (DTM), and Tc1/Mariner (DTT). 
They differ in their TIR length, and TSD length.


## How TIRs were Identified
Using exemplars from MTEC, ```target_mTEA_array.cacta.sh``` for example runs mTEA on each exemplar. 
Note that paths to perl libraries and mTEA should be updated, and that each script uses different TSD and TIR characteristics in the search.

## Assign Family Names

- ```cat *.tab > all_tir.tab``` to generate a file with all `*.tab` files

- ```Rscript tir_famnames.R``` uses `all_tir.tab` to reduce overlapping copies, and assign family names

