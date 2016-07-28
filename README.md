# maize_v4_TE_annotation
scripts and intermediate files used to annotate TEs in Jiao et al. 2016

----------------

__LTR Retrotransposons__

Scripts in ```nested_ltr```

This follows the structure of nested_ltr, here a submodule.

__SINEs__

Scripts in ```sine/```

Software needed:

- [SINE-Finder](http://www.plantcell.org/content/23/9/3117.full), [download](http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt) (This is a supplemental file at The Plant Cell; need to make executable, and rename to sine_finder.py)

	- I cannot make SINE-Finder function on reverse sequences. So I'm reporting SINEs only on the forward stand here, and will pick up sequences on the reverse strand with RepeatMasker.


__LINEs__

Scripts in ```line/```

Software needed:

- [MGEScanNonLTR](http://darwin.informatics.indiana.edu/cgi-bin/evolution/nonltr/nonltr.pl) I use the version generated for Galaxy [here](https://github.com/MGEScan/mgescan).

__TIR including MITEs__

Scripts in ```tir/```

Software needed:

- [MiteHunter](http://target.iplantcollaborative.org/mite_hunter.html), [download](http://target.iplantcollaborative.org/mite_hunter/MITE%20Hunter-11-2011.zip)

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

