import sys
import textwrap

def readFasta(filename):
	fastafile=open(filename, 'r')
	fastadict={}
	for line in fastafile:		
		if line.startswith('>') and line.strip().split(':')[0][1:] not in fastadict:
			seqname=line.strip().split(':')[0][1:]
			fastadict[seqname]=[]
		elif line.startswith('>') and line.strip().split(':')[0][1:] in fastadict:
#			seqname=line.strip().split(':')[0][1:]
			continue
		else:
			fastadict[seqname].append(line.strip())
	for entry in fastadict:
		fastadict[entry]=''.join(fastadict[entry])
	return fastadict

def printFasta(sequence, width=70):
	return '\n'.join(sequence[i:i+width] for i in range(0,len(sequence), width))

b=readFasta(sys.argv[1])
#out=open('stupid_test.txt', 'w')
#out.write('read in all the fasta entries')
for entry in b:
	print '>'+entry+'\n'+printFasta(''.join(b[entry]))


#contig=''
#seq=[]
#for line in fa:
#	if line.startswith('>'):
#		if contig == line.strip().split(':')[0][1:]:
#			continue
#		else:
#			if seq==[]:
#				print '>'+line.strip().split(':')[0][1:]
#				contig=line.strip().split(':')[0][1:]
#			else:
#				print textwrap.fill(''.join(seq), 70)
#				contig=line.strip().split(':')[0][1:]
#				print '>'+contig
##			print '>'+contig
#				seq=[]
#	else:
##		print line.strip()
#		seq.append(line.strip())
#
#print textwrap.fill(''.join(seq), 70)
