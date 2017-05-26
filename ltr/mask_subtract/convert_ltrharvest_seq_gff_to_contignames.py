import sys
import re
### takes gff from genometools, writes gff with contig names to stdout

# genome tools renames sequences so it can search by name
##     ##sequence-region   seq0 1 106474
##     #utg10007
##     seq0    LTRharvest      LTR_retrotransposon     58963   74561   .       ?       .       ID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=96.31;seq_number=0

gff=open(sys.argv[1], 'r')

seq=[]
contig=[]

linecount=0
for line in gff:
	fields=line.strip().split()
	if line.startswith('##sequence-region'):
		seq.append(line.split()[1])
	elif line.startswith('#') and re.search('#\d', line) is not None:
#	elif line.startswith('#B73'):   ### this is not general!
		contig.append(line.strip().split()[0][1:])
	elif len(fields)==9:
		contigname=contig[seq.index(fields[0])]
		fields[0]=contigname
		print '\t'.join(fields)
	else:
		print line.strip()

gff.close()
