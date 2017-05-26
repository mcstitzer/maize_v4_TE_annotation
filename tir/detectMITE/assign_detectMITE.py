import sys


#-----------------------------------------------------------
#>utg105632|820|967|11|4
#ACAGTCCGGTGCGGTGTCCGGTGCACCACTAAAATTCAACTCCGAATGCTGCGCTCTCGGGTTTCTGCGAAGGGAAAACACTGCCGCGGACCAGCCTGGCCCCACCTGGCAGATGGTGCACCGGACAGTCCGGTGCACACCGGACAGT
#>utg105667|224176|224323|11|4
#ACTGTCCGGTGTGCACCGGACTGTCCGGTGCACCCTCTGCCAGGTGGGGCCAGGCTGGTCCACGGTAGTGTTTTCCCTTCGCAGAAACCCGAGAGCGCAGCCTTGAGAGTTGAATTTTAGTGGCACACCGGACACCGCACCGGACTGT
#>utg130958|379626|379773|11|7
#ACTGTCCGGTGTGCACCGGACTGTCCGGTGCACCATCTACCAGGTGGGCCAGGCTGGCCCAGGGGAAGGCCTCCCCGGGCAGCAAAACCCGAGCGCGCCTGTTTCAGAAGTGAATTATAGTGGCGCACCGGACAGTGCACCAGACAGT
#>utg146178|1163996|1164143|11|4
#ACAGTCCGGTGCACTGTCCGGTGCGCCACTAAAATTCATTTCTGAAACTAGCGCTCTCGGGTTTCTTCGGGAGAAAACTCTTCCCCGGGGCCACCTTGGCCCCATCTAGTAGAGGGTGCACCGGACAGTCCGGTGCACACCGGACAGT
#


mitefa=open(sys.argv[1], 'r')

gffout=open(sys.argv[2], 'w')

out=0
tsd3=0
for line in mitefa:
	line=line.strip()
	if line.startswith('>'):
		out=0
		fields=line.strip().split('|')
		strand='*'  ## cannot know because these are mites without coding capacity
#####
		tsdlen=fields[4]
		tirlen=fields[3]
		if tsdlen == '2':
			print '>DTT_'+line[1:]
			out=1
			gffout.write(fields[0][1:]+'\tdetectMITE\tterminal_inverted_repeat_element\t'+fields[1]+'\t'+fields[2]+'\t.\t'+strand+'\t.\tID=DTT|'+line.strip()[1:]+';Name=TSDlength'+tsdlen+',TIRlength'+tirlen+'\n')
		elif tsdlen == '8':
			print '>DTA_'+line[1:]
			out=1
			gffout.write(fields[0][1:]+'\tdetectMITE\tterminal_inverted_repeat_element\t'+fields[1]+'\t'+fields[2]+'\t.\t'+strand+'\t.\tID=DTA|'+line.strip()[1:]+';Name=TSDlength'+tsdlen+',TIRlength'+tirlen+'\n')
		elif tsdlen == '3':
			cacta='>DTC_'+line[1:]
			pifharb='>DTH_'+line[1:]
			namefield=line.strip()[1:]
			tsd3=1
			out=1
		elif tsdlen =='9' or tsdlen=='10' or tsdlen=='11':
			print '>DTM_'+line[1:]
			out=1
			gffout.write(fields[0][1:]+'\tdetectMITE\tterminal_inverted_repeat_element\t'+fields[1]+'\t'+fields[2]+'\t.\t'+strand+'\t.\tID=DTM|'+line.strip()[1:]+';Name=TSDlength'+tsdlen+',TIRlength'+tirlen+'\n')
		else:
			print '>DTX_'+line[1:]
			out=1
			gffout.write(fields[0][1:]+'\tdetectMITE\tterminal_inverted_repeat_element\t'+fields[1]+'\t'+fields[2]+'\t.\t'+strand+'\t.\tID=DTX|'+line.strip()[1:]+';Name=TSDlength'+tsdlen+',TIRlength'+tirlen+'\n')
#			print line
			

#	else:
#		out=0
	elif line.startswith('#') or line.startswith('-'):
		continue
	else:
#		print str(out) + '_' + str(tsd3)
		if out==1 and tsd3==0: ### we've already printed the fasta header
			print line
		elif out==1 and tsd3==1:    ### we need to decide if this is cacta or pif harbinger
			if line.startswith('CACT'):
				print cacta
				gffout.write(fields[0][1:]+'\tdetectMITE\tterminal_inverted_repeat_element\t'+fields[1]+'\t'+fields[2]+'\t.\t'+strand+'\t.\tID=DTC|'+namefield+';Name=TSDlength'+tsdlen+',TIRlength'+tirlen+'\n')
			else:
				print pifharb
				gffout.write(fields[0][1:]+'\tdetectMITE\tterminal_inverted_repeat_element\t'+fields[1]+'\t'+fields[2]+'\t.\t'+strand+'\t.\tID=DTH|'+namefield+';Name=TSDlength'+tsdlen+',TIRlength'+tirlen+'\n')
			print line
			tsd3=0
		elif out==0:
			continue



