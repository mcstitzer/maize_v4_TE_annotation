I have run [detectMITE](https://www.nature.com/articles/srep19688) ([download](https://sourceforge.net/projects/detectmite/files/)) elsewhere because setting up the matlab dependencies is hard.

1. ```assign_detectMITE.py``` to put Wicker et al. (2007) superfamily designations on each TE, based on TSD length and TIR sequence (only for distinguishing CACTA and PIF/Harbinger which each have a 3bp TSD). 

	- ```python assign_detectMITE.py B73V4.both_pseudo_AND_unplaced.miteSet.fasta B73V4.both_pseudo_AND_unplaced.miteSet.sup.gff3 > B73V4.both_pseudo_AND_unplaced.miteSet.sup.fa```
	- this outputs a gff3 which will be imported in the next step, and outputs a fasta with superfamilies noted

2. ```cluster_families_MITE.sh``` to assign families based on Wicker et al. (2007)'s 808080 rule, and ```cluster_mtec.sh``` to associate MTEC family names with my family names.

3. ```generate_gff_assign_families_MITE.R``` will output a gff3 which assigns superfamily to each family, and gives TSD and TIR length in the Name section of field 9 of the gff3, and outputs a text file to allow relating fasta header to gff name (of the form DT?#####B73v4000001, where DT? is the superfamily three character designation, ##### is the family number, and B73v4000001 is the sequential number of this copy). 


__NOTE__

DetectMITE outputs some overlapping TEs. I'm reporting all, which may be a bad thing if you're wanting to identify genomic coverage or copy number. 
