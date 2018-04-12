import sys

### how to run
#srun -p bigmemh python helitron_scanner_out_to_tabout.py B73.Mhap2.quiver.fasta.HelitronScanner.draw.hel.fa > B73.Mhap2.quiver.fasta.HelitronScanner.gff3
#srun -p bigmemh python helitron_scanner_out_to_tabout.py B73.Mhap2.quiver.fasta.HelitronScanner.draw.rc.hel.fa >> B73.Mhap2.quiver.fasta.HelitronScanner.gff3




hsfa=open(sys.argv[1], 'r')

## be careful if you're rerunning because this appends
outfile=open(sys.argv[2], 'a')

### example header line
#>10_#SUB_75503-83673 [forward] 8171 bps; scores=10:6 Multi_5'_ends:
#>10_#SUB_3564268-3578686 [forward] 14419 bps; scores=8:7 Multi_5'_ends:3573243:8

## want to keep track of the first and last couple of bp to get a crude estimate of family
## and to predict hairpins in the 3' end

seq='first'
for line in hsfa:
	if line.startswith('>'):
		fields=line.strip().split()
		posn=fields[0].split('_')
		### because fields can't be split by underscore easily, as contig names have these, generate an offset
		try:
			offset=posn.index('#SUB')
		except ValueError:
			offset=0
		chrom='_'.join(posn[0:offset])[1:]
#		start=posn[2].split('-')[0]
#		end=posn[2].split('-')[1]
		orientation=fields[1][1:-1]
		if orientation=='forward':
			orientation='+'
			start=posn[1+offset].split('-')[0]
			end=posn[1+offset].split('-')[1]
		elif orientation=='reverse':
			orientation='-'
			end=posn[1+offset].split('-')[0]  ## need to switch because gff always expects start<end
			start=posn[1+offset].split('-')[1]

		length=fields[2]
		scores=fields[4][7:]
		multiend=fields[5][14:]
		if seq!='':
#			print chrom+'\t'+start+'\t'+end+'\t'+orientation+'\t'+length+'\t'+scores.split(':')[0]+'\t'+scores.split(':')[1]+'\t'+multiend+'\t'+seq[:30]+'\t'+seq[-30:]
#			print chrom+'\tHelitronScanner\thelitron\t'+start+'\t'+end+'\t.\t'+orientation+'\t.\t'+'ID='+chrom+'_'+start+'_'+end+';Name=LeftScore'+scores.split(':')[0]+',RightScore'+scores.split(':')[1]
			print chrom+'\t'+start+'\t'+end+'\t'+orientation+'\t'+chrom+'_'+start+'_'+end+'\t'+scores.split(':')[0]+'\t'+scores.split(':')[1]+'\t'+seq[:30]+'\t'+seq[-30:]
			seq=''
		outfile.write('>'+chrom+'_'+start+'_'+end+'\n')
	else:
		seq+=line.strip()  ### i know this is bad because it generates a new string each time.
		outfile.write(line)

## to get the last one
#print chrom+'\t'+start+'\t'+end+'\t'+orientation+'\t'+length+'\t'+scores.split(':')[0]+'\t'+scores.split(':')[1]+'\t'+multiend+'\t'+seq[:30]+'\t'+seq[-30:]
### gff format
#print chrom+'\tHelitronScanner\thelitron\t'+start+'\t'+end+'\t.\t'+orientation+'\t.\t'+'ID='+chrom+'_'+start+'_'+end+';Name=LeftScore'+scores.split(':')[0]+',RightScore'+scores.split(':')[1]
##### tabout
print chrom+'\t'+start+'\t'+end+'\t'+orientation+'\t'+chrom+'_'+start+'_'+end+'\t'+scores.split(':')[0]+'\t'+scores.split(':')[1]+'\t'+seq[:30]+'\t'+seq[-30:]
