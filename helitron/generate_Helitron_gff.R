library(rtracklayer)
library(stringr)
library(plyr)


GENOMENAME=commandArgs(trailingOnly=TRUE)
#GENOMENAME='B73V4.both_pseudo_AND_unplaced'

######### Read in silix results
a=read.table(paste('families/', GENOMENAME, '.fa.HelitronScanner.draw.bothorientations.tabnames.final30bp.8080.fnodes', sep=''))
a$V1=as.character(a$V1)
a$V2=as.character(a$V2)


######### Read in helitronscanner results tabout, fix naming, classes
hel=read.table(paste(GENOMENAME, '.fa.HelitronScanner.draw.bothorientations.tabout', sep=''), header=T, sep='\t')

### example line 10      299043  314401  +       10_299043_314401        7       7       TCTATATATACATATTACTCCATGTCTATA  AAGCTTCTCTACCATGGTGGTGCTTGCTAA
names(hel)=c('chr', 'start', 'end', 'orientation', 'ID', 'score5', 'score3', 'fivebp', 'threebp')

hel$family=mapvalues(hel$ID, from=a$V2, to=a$V1)

### assign mtec families

pm=read.table(paste('families/', GENOMENAME, '.TEDB.8080.searchglobal.toponly.out', sep=''), header=F)
pm$wicker=str_split_fixed(pm$V1, '_', 3)[,1]
pm=pm[pm$wicker=='DHH',]
pm$family=str_split_fixed(pm$V1, '_', 3)[,2]
pm.gd=merge(pm, a, by.x='V2', by.y='V2', all.x=T)

#get an mtec family if any of the family members of mine matched
hel$mtec.fam=NA
for (i in 1:nrow(pm.gd){
	hel$mtec.fam[hel$family==pm.gd$V1.y[i]]=pm.gd$family[i]
	}

### these are all DHH in wicker scheme, so assign all superfamily DHH
### sort by biggest family first
rankedfams=names(rev(sort(table(hel$family))))
hel$famrank=match(hel$family, rankedfams)
hel$rankedfamname=sapply(1:nrow(hel), function(x) paste('HEL0', str_pad(hel$famrank[x], 5, pad='0'), sep=''))



hel.gr=GRanges(seqnames=hel$chr, ranges=IRanges(start=hel$start, end=hel$end), strand=hel$orientation)
mcols(hel.gr)$ID=hel$ID
mcols(hel.gr)$score5=hel$score5
mcols(hel.gr)$score3=hel$score3
mcols(hel.gr)$rankedfamname=hel$rankedfamname
mcols(hel.gr)$mtec.fam=hel$mtec.fam

f=merge(hel, a, by.x='ID', by.y='V2', all=T)
f$fam=table(f$V1)[f$V1]
mcols(hel.gr)$fam=f$fam


## output gff

