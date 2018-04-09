# maize_v4_TE_annotation
scripts and intermediate files used to annotate TEs in [Jiao et al. 2016](http://biorxiv.org/content/early/2017/04/20/079004)

----------------

__LTR Retrotransposons__

scripts in ```ltr```

Software needed:

- ncbi blast+

- [genometools](http://www.genometools.org), ([download](http://www.genometools.org/pub/genometools-1.5.7.tar.gz)), need to pass `64bit=yes with-hmer=yes threads=yes` to make, make install for ltrdigest hmm searches in parallel. I also had to pass `cairo=no` as well because I didn't have the right cairo libraries and it wouldn't compile otherwise

- [silix](http://lbbe.univ-lyon1.fr/Download,3009.html?lang=en), ([download](ftp://pbil.univ-lyon1.fr/pub/logiciel/silix//silix-1.2.10.tar.gz)), need to compile with ```--enable-mpi``` and ```--enable-verbose```

- hmmer (genometools with download and compile hmmer2 if you run ```make with-hmmer=yes```)

Files needed, can be downloaded by `get_tRNA_hmm_dbs.sh` in ```ltr``` directory:

- download [hmms](http://gydb.org/gydbModules/collection/collection/db/GyDB_collection.zip) of TE protein coding domains from gydb in directory `gydb_hmms`, will be used to identify protein coding domains of TE models
	
	-need to fix a hmm with name ty1/copia because this is used as a filename by ltrdigest. to remove the forward slash:  ```sed -i "s#ty1/copia#ty1-copia#g" gydb_hmms/GyDB_collection/profiles/AP_ty1copia.hmm```

- download [tRNAs](http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz) of all eukaryotes

__SINEs__

Scripts in ```sine/```

Software needed:

- [SINE-Finder](http://www.plantcell.org/content/23/9/3117.full), [download](http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt) (This is a supplemental file at The Plant Cell; need to make executable, and rename to sine_finder.py)

	- I cannot make SINE-Finder function on reverse sequences. So I'm reporting SINEs only on the forward stand here, and will pick up sequences on the reverse strand with RepeatMasker.


__LINEs__

Scripts in ```line/```

Software needed:

- [MGEScanNonLTR](http://darwin.informatics.indiana.edu/cgi-bin/evolution/nonltr/nonltr.pl) I use the version generated for Galaxy [here](https://github.com/MGEScan/mgescan).

- [TARGeT](https://academic.oup.com/nar/article/37/11/e78/1094076), [download](https://github.com/BradCavinder/TARGeT), also used for TIRs


__TIR including MITEs__

Scripts in ```tir/```

Software needed:

- [detectMITE](https://www.nature.com/articles/srep19688), [download](https://sourceforge.net/projects/detectmite/)

- [TARGeT](https://academic.oup.com/nar/article/37/11/e78/1094076), [download](https://github.com/BradCavinder/TARGeT)

- [mTEA](https://github.com/stajichlab/mTEA), genometools (see above, already installed for ltr annotation) 

  - mTEA needs [fasta36](http://faculty.virginia.edu/wrpearson/fasta/fasta36/) (specifically ggsearch36), bioperl, blast, muscle, supplied blogo directories to be put into PERL5LIB and PATH


__Helitrons__

Scripts in ```helitron/```

Software needed:

- [HelitronScanner](http://bo.csam.montclair.edu/du/software/helitronscanner)

------------


**Finding Homologous Fragments from Degraded TEs**
----------------

Software needed:

- [RepeatMasker](http://www.repeatmasker.org/), with prerequisites [here](http://www.repeatmasker.org/RMDownload.html)

