library(rtracklayer)
library(stringr)


#GENOMENAME='B73V4.both_pseudo_AND_unplaced'
GENOMENAME=commandArgs(trailingOnly=TRUE)

######### Read in silix results (family assignments)
a=read.table(paste(GENOMENAME, '-matches.noTSD.8080.fnodes', sep=''))

### Switch SINE to RST, the wicker 3 letter code
## because these are plants, and found with the rna pol A and B boxes, these are all RST
## remove one digit because there are not 100,000 families, one order of mag less
#f$V1=sub('SINE0', 'RST', f$V1)  ### let's switch and do as early as we can!!!
a$V1=sub('SINE0', 'RST', a$V1)


######### Read in sine_finder results, fix naming, classes
sine=read.table(paste(GENOMENAME, '-matches.csv', sep=''), header=T)
sine$sequence=NULL
sine$direct='+'
sine$namecomp=paste(sine$name, sine$start, sine$end, 'F', paste('TSDlen', sine$TSD.len, sep=''), paste('TSDscore', sine$TSD.score, sep=''), paste('TSDmism', sine$TSD.mism, sep=''), sep='_')

f=merge(sine, a, by.x='namecomp', by.y='V2', all=T)
f$fam=table(f$V1)[f$V1] 


### give each TE copy a unique identifier
f$Name=''
### can't figure out how to assign properly, so giving up and doing a for loop
#f$Name=sapply(names(table(f$V1)), function(x) f$Name[f$V1==x]=paste(x, 'B73v4', str_pad(1:sum(f$V1==x), 5, pad='0'), sep=''))
for (x in names(table(f$V1))){
	f$Name[f$V1==x]=paste(x, 'B73v4', str_pad(1:sum(f$V1==x), 5, pad='0'), sep='')
	}
	


### match to entries in the TE database (maizetedb.org)
pm=read.table(paste(GENOMENAME, '-matches.noTSD.TEDB8080.out', sep=''), header=F)
pm$wicker3=str_split_fixed(pm$V1, '_', 3)[,1]
pm=pm[pm$wicker3=='RST',]
#pm.gd=f$namecomp %in% pm$V2
pm$family=str_split_fixed(pm$V1, '_', 3)[,2]

pm.gd=merge(pm, a, by.x='V2', by.y='V2', all.x=T)


f$mtec.fam=NA
for (i in 1:nrow(pm.gd)){
	f$mtec.fam[f$V1==pm.gd$V1.y[i]]=pm.gd$family[i]
	}

## output the gff

f.gff=data.frame(f$name, 'SineFinder', 'SINE_element', f$start, f$end, '.', '+', '.', paste('ID=', f$Name, ';Name=', paste(f$Name, f$mtec.fam, paste('TSDlen', f$TSD.len, sep=''), paste('TSDmismat', f$TSD.mism, sep=''), sep='_'), sep=''))
write.table(f.gff, paste(GENOMENAME, '.RST.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

### also keep track of fasta names and newly assigned gffnames as to easily convert between the two (e.g. switching fasta names)
write.table(data.frame(f$Name, f$namecomp), paste(GENOMENAME, '.RST.gffname.fastaname.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

