import sys

infile=open(sys.argv[1], 'r')

tes=[]
for line in infile:
	tes.append(line.strip())

infile.close()

fafile=open(sys.argv[2], 'r')
for line in fafile:
	line=line.strip()
	if line.startswith('>'):    ### check if it is a te we want
		if line[1:] in tes:
			print line
			useTE=True
		else:
			useTE=False
	else:
		if useTE==True:
			print line
		else:
			continue

fafile.close()
