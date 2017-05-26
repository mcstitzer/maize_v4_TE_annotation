library(rtracklayer)
library(stringr)
library(plyr)

GENOMENAME='B73V4.both_pseudo_AND_unplaced'

######### Read in silix results
a=read.table(paste(GENOMENAME, '.miteSet.8080.fnodes', sep=''))
a$V1=as.character(a$V1)
a$V2=as.character(a$V2)

#a$start=str_split_fixed(a$V2, '\\|', 6)[,2]
#a$end=str_split_fixed(a$V2, '\\|', 6)[,3]
#a$sequence=str_split_fixed(a$V2, '\\|', 6)[,1]



######### Read in detectMITE results gff, fix naming, classes
mite=import.gff(paste(GENOMENAME, '.miteSet.sup.gff3', sep=''))
mcols(mite)$sup=str_split_fixed(mcols(mite)$ID, '\\|', 6)[,1]
mcols(mite)$V2=substring(mcols(mite)$ID, 5)
mcols(mite)$TSDlen=str_split_fixed(mcols(mite)$ID, '\\|', 6)[,6]
mcols(mite)$TIRlen=str_split_fixed(mcols(mite)$ID, '\\|', 6)[,5]
mcols(mite)$family=mapvalues(mcols(mite)$V2, from=a$V2, to=a$V1)
mcols(mite)$famsize=table(mcols(mite)$family)[mcols(mite)$family]

	

pm=read.table(paste(GENOMENAME, '.miteSet.noSep.TEDB.8080.searchglobal.toponly.out', sep=''), header=F)
pm$wicker3=str_split_fixed(pm$V1, '_', 3)[,1]
pm=pm[pm$wicker3 %in% c('DTT', 'DTA', 'DTH', 'DTM', 'DTC', 'DTX'),]
pm$family=str_split_fixed(pm$V1, '_', 3)[,2]

pm.gd=merge(pm, a, by.x='V2', by.y='V2', all.x=T)

mcols(mite)$mtec.sup=NA
mcols(mite)$mtec.fam=NA
mcols(mite)$mtec=NA
for (i in 1:nrow(pm.gd)){
	mcols(mite)$mtec.sup[mcols(mite)$family==pm.gd$V1.y[i]]=pm.gd$wicker3[i]
	mcols(mite)$mtec.fam[mcols(mite)$family==pm.gd$V1.y[i]]=pm.gd$family[i]
	mcols(mite)$mtec[mcols(mite)$family==pm.gd$V1.y[i]]=paste(pm.gd$wicker3[i], pm.gd$family[i], sep='')
	}


## assign superfamily to each family based on whether protein order can be assigned
fam.sup=sapply(names(table(mcols(mite)$family)), function(x) names(which.max(table(mcols(mite)$sup[which(mcols(mite)$family==x)]))))
mcols(mite)$fam.sup=mapvalues(mcols(mite)$family, from=names(unlist(fam.sup)), to=unlist(fam.sup))

## sometimes family names are poorly called based on 
## this mostly rescues some large families, and occurs when the TSD or TIR is a bit off (could happen with insertion preference, like it's much easier for a DTT two base pair TSD to be mistaken as a longer TSD by chance)
mite$rescue=mite$fam.sup
mite$rescue[mite$fam.sup=='DTX']=mite$mtec.sup[mite$fam.sup=='DTX']
mite$rescue[is.na(mite$rescue)]='DTX'

### switch the naming of families from silix to one that reflects superfamily ala Wicker et al. (2007), and is 5 digits long
mcols(mite)$famname=sapply(1:length(mite), function(x) sub('MITE0', mcols(mite)$rescue[x], mcols(mite)$family[x]))

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
mcols(mite)$Name=NULL
for (x in names(table(mcols(mite)$famname))){
  mcols(mite)$Name[mcols(mite)$famname==x & !is.na(mcols(mite)$famname)]=paste(x, 'B73v4', str_pad(1:sum(mcols(mite)$famname[!is.na(mcols(mite)$famname)]==x), 5, pad='0'), sep='')
}



## now I can output the gff!

mite.gff=data.frame(seqnames(mite), 'detectMITE', 'terminal_inverted_repeat_element', start(mite), end(mite), '.', strand(mite), '.', paste('ID=', mcols(mite)$Name, ';Name=', paste(mcols(mite)$Name, mcols(mite)$mtec, paste('TSDlen', mcols(mite)$TSDlen, sep=''), paste('TIRlen', mcols(mite)$TIRlen, sep=''), sep='_'), sep=''))

write.table(mite.gff, paste(GENOMENAME, '.detectMITE.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)


### also keep track of fasta names and newly assigned gffnames as to easily convert between the two (e.g. switching fasta names)
write.table(data.frame(mcols(mite)$Name, mcols(mite)$V2, mcols(mite)$ID), paste(GENOMENAME, '.detectMITE.gffname.fastaname.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

