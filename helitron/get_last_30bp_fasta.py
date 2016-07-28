import sys


def readFasta(filename):
	fastafile=open(filename, 'r')
	fastadict={}
	for line in fastafile:		
		if line.startswith('>'):
			seqname=line.strip()[1:]
			fastadict[seqname]=[]
		else:
			fastadict[seqname].append(line.strip())
	for entry in fastadict:
		fastadict[entry]=''.join(fastadict[entry])
	return fastadict


fdict=readFasta(sys.argv[1])

for entry in fdict:
	print '>'+entry+'\n'+fdict[entry][-30:]


