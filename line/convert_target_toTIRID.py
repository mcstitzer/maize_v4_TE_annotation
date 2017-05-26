
import sys

fa=open(sys.argv[1], 'r')

for line in fa:
	if line.startswith('>'):
		fields=line.strip().split()
#		print fields
		print '>'+fields[1].split(':')[1]+'_'+fields[2].split(':')[1]+':'+fields[4].split('(')[1]+'..'+fields[6][:-1]
	else:
		print line.strip()

