

import sys

infile=open(sys.argv[1], 'r')
outfile=open(sys.argv[2],'w')


#>10 F 1960983:1961124 TSD-len=12;TSD-score=11;TSD-mism=1
#TAAAATTTAAATgggccagGCTGGctcgaaatttaaacggatcgggcttttttgggcttgagccgagccggacgGTTCGAatatacacctataatctaaaatatggggaaaaagcttcgttgtAAAAAACAAAATATAAAT

fastadict={}

for line in infile:
	### READ FASTA
	line=line.strip()
	if line.startswith('\n'):
		continue
	if line.startswith('>'):
		dictkey=line[1:]
		fastadict[dictkey]=[]
	else:
		fastadict[dictkey].append(line)
for key in fastadict:
	fastadict[key]=''.join(fastadict[key])


### remove the tsd from the fasta entry
### should eventually store the removed tsd sequences in the fasta header


### I'm going to also add in some underscores to this disaster of a fasta header
#>10 F 1960983:1961124 TSD-len=12;TSD-score=11;TSD-mism=1

notsddict={}

for key in fastadict:
	tsdlen=key.split('=')[1]
	tsdlen=tsdlen.split(';')[0]
	tsdlen=int(tsdlen)
	notsddict[key]=fastadict[key][tsdlen:len(fastadict[key])-tsdlen]
	newname=key.split()
#	outfile.write('>'+key+'\n'+notsddict[key]+'\n')
	outfile.write('>'+newname[0]+'_'+newname[2].split(':')[0]+'_'+newname[2].split(':')[1]+'_'+newname[1]+'_TSDlen'+newname[3].split('=')[1].split(';')[0]+'_TSDscore'+newname[3].split('=')[2].split(';')[0]+'_TSDmism'+newname[3].split('=')[3]+'\n'+notsddict[key]+'\n')

