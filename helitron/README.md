## Helitron background

Helitrons are loosely structurally defined, as their 5' and 3' terminii are a handful of nucleotides long. 
We apply HelitronScanner, but caution that many of these elements are candidates, and should not be believed to be true helitron elements without further study.

## How Helitrons Were Identified

The script ```run_helitron_scanner.sh``` identifies direct and reverse complement orientation helitrons, converts the HelitronScanner output to a tabular format, grabs the final 30 bp of the element for clustering, does the clustering, assigns families, and outputs a gff.

Note that to read in entire maize chromosomes (and keep coordinates safe), HelitronScanner requires a lot of memory (I gave it 300 Gb, but I may have been overdoing it).
