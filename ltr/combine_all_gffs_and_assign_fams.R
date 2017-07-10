library(rtracklayer)
library(stringr)
library(plyr)


GENOMENAME='B73V4.both_pseudo_AND_unplaced'
MAXSUBTRACT=80
MAXUNPLACED=9
UNPLACEDGENOME='B73V4.unplaced'
PSEUDOMOLECULE='B73V4.pseudomolecule'

#### GET LTRS.MASK AND PROT.MASK

### get ltrharvest results - need to do chromosomes and unplaced separately because i did them separately (and won't again!)

masked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '/', PSEUDOMOLECULE, '.subtract', 1:MAXSUBTRACT, '.ltrharvest.contignames.gff3.contigpositions.gff3', sep='')
filenames=c(paste('mask_subtract/',PSEUDOMOLECULE, '.ltrharvest.contignames.gff3', sep=''), masked)


ltrs=lapply(filenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
                                    mcols(a)$nest=which(filenames==x)-1
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    mcols(a)$ltr_similarity=as.numeric(mcols(a)$ltr_similarity)
                                    mcols(a)$nestlevel=which(filenames==x)-1
                                    return (a)
                                    }
            )
            
ltrs.mask=do.call(c, ltrs)

## do for unplaced ones
unplaced=paste('mask_subtract_unplaced/subtract',1:MAXUNPLACED, '/', UNPLACEDGENOME, '.subtract', 1:MAXUNPLACED, '.ltrharvest.sorted.contignames.gff3.contigpositions.gff3', sep='')
unplacedfilenames=c(paste('mask_subtract_unplaced/', UNPLACEDGENOME, '.ltrharvest.sorted.contignames.gff3', sep=''), unplaced)

unplacedltrs=lapply(unplacedfilenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
                                    mcols(a)$nest=which(unplacedfilenames==x)-1
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    mcols(a)$ltr_similarity=as.numeric(mcols(a)$ltr_similarity)
                                    mcols(a)$nestlevel=which(unplacedfilenames==x)-1
                                    return (a)
                                    }
            )
ltrs.mask=c(ltrs.mask, do.call(c, unplacedltrs))  ### and we expect these to have no sequence levels in common, so the warning is fine

### add in the coverage of each fragment - how many are overlapping
seqcov=coverage(ltrs.mask)
mcols(ltrs.mask)$cov=NA  ## initialize the variable before loop
for (i in as.character(seqlevels(ltrs.mask))){ ltrs.mask[as.character(seqnames(ltrs.mask))==i,]$cov=as.numeric(seqcov[[i]][start(ltrs.mask[as.character(seqnames(ltrs.mask))==i,])])}



#### get ltrdigest tabout files (special here because i ran the unplaced separately from the chromosomes (should not have!)
##### MAKE SURE THESE HAVE THE ^B REMOVED!

dignames=paste('mask_subtract/digest_tabout_vim/', PSEUDOMOLECULE, '.subtract', 1:MAXSUBTRACT, '.ltrdigest_tabout.csv', sep='')
digfiles=c(paste('ltrdigest-redo/B73V4.pseudomolecule.ltrdigest_tabout.csv', sep=''), dignames)


split.fam=function(s) strsplit(as.character(s), '_')[[1]][1]
prot=lapply(digfiles, function(x) {m=read.table(x, sep='\t', header=T)
                                   m$sequence=str_split_fixed(as.character(m$sequence), '_', 3)[,1]
                                   m$rle=''
                                   print(x)
                                   m$rle[m$Protein.domain.hits!='']=sapply(which(m$Protein.domain.hits!=''), function(x)  rle(sapply(ldply(strsplit(as.character(m$Protein.domain.hits[x]), '/')), split.fam))$values) 
                                   m$geneorder=sapply(1:nrow(m), function(x) paste(m$rle[[x]], collapse='.'))
                                   m$nestlevel=which(digfiles==x)-1
                                   return(m)
                                   }
            )
prot.mask=do.call(rbind, prot)

##### unplaced now
## same order as harvests
digunplaced=paste('mask_subtract_unplaced/digest_tabout_vim/', UNPLACEDGENOME, '.subtract', 1:MAXUNPLACED, '.ltrdigest_tabout.csv', sep='')
unplaceddigfiles=c(paste('ltrdigest_unplaced/', UNPLACEDGENOME, '.ltrdigest_tabout.csv', sep=''), digunplaced)
unplacedprot=lapply(unplaceddigfiles, function(x) {m=read.table(x, sep='\t', header=T)
                                   m$sequence=paste(str_split_fixed(as.character(m$sequence), '_', 3)[,1], str_split_fixed(as.character(m$sequence), '_', 3)[,2], sep='_')
                                   m$rle=''
                                   print(x)
                                   m$rle[m$Protein.domain.hits!='']=sapply(which(m$Protein.domain.hits!=''), function(x)  rle(sapply(ldply(strsplit(as.character(m$Protein.domain.hits[x]), '/')), split.fam))$values) 
                                   m$geneorder=sapply(1:nrow(m), function(x) paste(m$rle[[x]], collapse='.'))
                                   m$nestlevel=which(unplaceddigfiles==x)-1
                                   return(m)
                                   }
            )
prot.mask=rbind(prot.mask, do.call(rbind, unplacedprot))



############# get the same for the ltrdigest gff3's so that I can add strand
dmasked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '/', PSEUDOMOLECULE, '.subtract', 1:MAXSUBTRACT, '.ltrdigest.gff3', sep='')
dfilenames=c(paste(PSEUDOMOLECULE, '.ltrdigest-redo.gff3', sep=''), dmasked)

dltrs=lapply(dfilenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
									print(x)
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    ltr_similarity=as.numeric(mcols(a)$ltr_similarity) ## can use this to confirm same order
                                    nestlevel=which(dfilenames==x)-1
                                    ID=mcols(a)$ID
                                    mcols(a)=NULL ### get rid of this because what comes first determines columns, which differ between file
                                    mcols(a)$ltr_similarity=ltr_similarity
                                    mcols(a)$nestlevel=nestlevel
                                    mcols(a)$ID=ID
                                    return (a)
                                    }
            )
            
dltrs.mask=do.call(c, dltrs)

## do for unplaced ones
dunplaced=paste('mask_subtract_unplaced/subtract',1:MAXUNPLACED, '/', UNPLACEDGENOME, '.subtract', 1:MAXUNPLACED, '.ltrdigest.gff3', sep='')
dunplacedfilenames=c(paste(UNPLACEDGENOME, '.ltrdigest.gff3', sep=''), dunplaced)

dunplacedltrs=lapply(dunplacedfilenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
									print(x)
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    ltr_similarity=as.numeric(mcols(a)$ltr_similarity) ## can use this to confirm same order
                                    nestlevel=which(dunplacedfilenames==x)-1
                                    ID=mcols(a)$ID
                                    mcols(a)=NULL ### get rid of this because what comes first determines columns, which differ between file
                                    mcols(a)$ltr_similarity=ltr_similarity
                                    mcols(a)$nestlevel=nestlevel
                                    mcols(a)$ID=ID
                                    return (a)
                                    }
            )
dltrs.mask=c(dltrs.mask, do.call(c, dunplacedltrs))  ### and we expect these to have no sequence levels in common, so the warning is fine



########## OH NO THESE ARE THE WRONG ORDER!!!!!######
###### these will be good to switch the order, as the start and end positions are the subtracted coordinates

### get ltrharvest results - need to do chromosomes and unplaced separately because i did them separately (and won't again!)

smasked=paste('mask_subtract/subtract', 1:MAXSUBTRACT, '/', PSEUDOMOLECULE, '.subtract', 1:MAXSUBTRACT, '.ltrharvest.gff3', sep='')
sfilenames=c(paste(PSEUDOMOLECULE, '.ltrharvest.gff3', sep=''), smasked)


sltrs=lapply(sfilenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
									print(x)
                                    mcols(a)$nest=which(sfilenames==x)-1
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    mcols(a)$ltr_similarity=as.numeric(mcols(a)$ltr_similarity)
                                    mcols(a)$nestlevel=which(sfilenames==x)-1
                                    return (a)
                                    }
            )
            
sltrs.mask=do.call(c, sltrs)

## do for unplaced ones
sunplaced=paste('mask_subtract_unplaced/subtract',1:MAXUNPLACED, '/', UNPLACEDGENOME, '.subtract', 1:MAXUNPLACED, '.ltrharvest.gff3', sep='')
sunplacedfilenames=c(paste(UNPLACEDGENOME, '.ltrharvest.gff3', sep=''), sunplaced)

sunplacedltrs=lapply(sunplacedfilenames, function(x) {a=import.gff(x, version='3', asRangedData=F)
                                    mcols(a)$nest=which(sunplacedfilenames==x)-1
                                    a=a[mcols(a)$type=='LTR_retrotransposon',]
                                    mcols(a)$ltr_similarity=as.numeric(mcols(a)$ltr_similarity)
                                    mcols(a)$nestlevel=which(sunplacedfilenames==x)-1
                                    return (a)
                                    }
            )
sltrs.mask=c(sltrs.mask, do.call(c, sunplacedltrs))  ### and we expect these to have no sequence levels in common, so the warning is fine

#########

## confirm these are in the same order

all(mcols(ltrs.mask)$ltr_similarity==mcols(dltrs.mask)$ltr_similarity)   ## concerning that this is not true
all(mcols(ltrs.mask)$ID==mcols(dltrs.mask)$ID)   ### but this has to be right! does this mean i'm switching contig names wrong?

##fine now
### there are different orders in the digest and harvest and i don't know why but i'm solving it.
#mcols(ltrs.mask)$compare=paste('seq', mcols(ltrs.mask)$seq_number, mcols(ltrs.mask)$ltr_similarity, mcols(ltrs.mask)$nestlevel, sep='')
#mcols(dltrs.mask)$compare=paste(seqnames(dltrs.mask), mcols(dltrs.mask)$ltr_similarity, mcols(dltrs.mask)$nestlevel, sep='')
#ltrs.mask=ltrs.mask[order(match(mcols(ltrs.mask)$compare, mcols(dltrs.mask)$compare)),]




## these should be the same order as the ltrs.mask object
## also check that the start values are the same to be super sure
length(ltrs.mask)==nrow(prot.mask)
all(start(ltrs.mask)==prot.mask$element.start) ### this will be false now that these are subtracted!!!!
all(as.character(seqnames(ltrs.mask))==prot.mask$sequence)
all(as.character(seqnames(ltrs.mask))==prot.mask$sequence)


## if these are all fine, it's time to assign strand to the ltrs!

strand(ltrs.mask)=strand(dltrs.mask)





##############################################
##  Read in and process family assignments  ##
##############################################

######### Read in silix results
a=read.table(paste('families/', GENOMENAME, '.ltrdigest_5ltr.combined.8080.fnodes', sep=''))
a$V1=as.character(a$V1)
a$V2=as.character(a$V2)

### get switch positions for chromosome to subtracted chromosome here
glen=readRDS(paste('mask_subtract/', PSEUDOMOLECULE, '_gl_', MAXSUBTRACT, '.RDS', sep=''))

### read in each subtracted family file, to then switch their positions
filenames=paste('families/subtracted/', PSEUDOMOLECULE, '.ltrdigest.subtract', 1:MAXSUBTRACT, '.ORIG5pLTR.id80.cov80.out', sep='')


### switch the positions in each file, write an updated file with .contigpositions.famnames.out as the suffix
for (sublev in 1:length(filenames)){
     b6=read.table(filenames[sublev], sep='\t', header=F)
     b6$V2=as.character(b6$V2)
     posns=str_split_fixed(b6$V1, '_', 3)
     posns=as.data.frame(posns)
     posns$V2=as.numeric(as.character(posns$V2))
     posns$V3=as.numeric(as.character(posns$V3))
     seqnames=names(table(posns[,1]))
     for (i in seqnames){
          rowindex=as.vector(posns[,1]==i)   ## only want to change positions in the gff for this chromosome/contig
          chrlist=glen[[i]]    ### these are the gff number at which each bp was annotated as TE
          chrlist.subtract=which(is.na(chrlist) | chrlist>sublev)   ### JUST GREATER THAN BECAUSE NOT STARTING WITH ORIGINAL HARVEST FILE these are the positions that would have been visible to ltrharvest when annotating this TE
          posns[rowindex,2]=chrlist.subtract[posns[rowindex,2]]   ### do end first so genomic ranges doesn't complain about negative widths
          posns[rowindex,3]=chrlist.subtract[posns[rowindex,3]]
          }
     b6$V1=paste(posns$V1, posns$V2, posns$V3, sep='_')
     b6$V2=mapvalues(b6$V2, from=a$V2, to=a$V1, warn_missing=F)
     write.table(b6, paste(filenames[sublev], '.contigpositions.famnames.out', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
     }


### REPEAT WITH UNPLACED CONTIGS!
### get switch positions for chromosome to subtracted chromosome here
glen=readRDS(paste('mask_subtract_unplaced/', UNPLACEDGENOME, '_gl_', MAXUNPLACED, '.RDS', sep=''))
unplacedfilenames=paste('families/subtracted/', UNPLACEDGENOME, '.ltrdigest.subtract', 1:MAXUNPLACED, '.ORIG5pLTR.id80.cov80.out', sep='')


### switch the positions in each file, write an updated file with .contigpositions.famnames.out as the suffix
for (sublev in 1:length(unplacedfilenames)){
     b6=read.table(unplacedfilenames[sublev], sep='\t', header=F)
     b6$V2=as.character(b6$V2)
     posns=str_split_fixed(b6$V1, '_', 4)
     posns=as.data.frame(posns)
     posns$V3=as.numeric(as.character(posns$V3))    ### this is different because contig names have an underscore!
     posns$V4=as.numeric(as.character(posns$V4))
     seqnames=names(table(paste(posns[,1], posns[,2], sep='_')))
     for (i in seqnames){
          rowindex=as.vector(posns[,1]==i)   ## only want to change positions in the gff for this chromosome/contig
          chrlist=glen[[i]]    ### these are the gff number at which each bp was annotated as TE
          chrlist.subtract=which(is.na(chrlist) | chrlist>sublev)   ### JUST GREATER THAN BECAUSE NOT STARTING WITH ORIGINAL HARVEST FILE these are the positions that would have been visible to ltrharvest when annotating this TE
          posns[rowindex,2]=chrlist.subtract[posns[rowindex,2]]   ### do end first so genomic ranges doesn't complain about negative widths
          posns[rowindex,3]=chrlist.subtract[posns[rowindex,3]]
          }
     b6$V1=paste(posns$V1, posns$V2, posns$V3, sep='_')
     b6$V2=mapvalues(b6$V2, from=a$V2, to=a$V1, warn_missing=F)
     write.table(b6, paste(unplacedfilenames[sublev], '.contigpositions.famnames.out', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
     }






### now read them back in, could fix this to be more efficient, i know
famfiles=paste('families/subtracted/', PSEUDOMOLECULE, '.ltrdigest.subtract', 1:MAXSUBTRACT, '.ORIG5pLTR.id80.cov80.out.contigpositions.famnames.out', sep='')
## also unplaced
famfiles=c(famfiles, paste('families/subtracted/', UNPLACEDGENOME, '.ltrdigest.subtract', 1:MAXUNPLACED, '.ORIG5pLTR.id80.cov80.out.contigpositions.famnames.out', sep=''))

fams=lapply(famfiles, function(x) {a=read.table(x, sep='\t', header=F)
									a$V1=as.character(a$V1)
									a$V2=as.character(a$V2)
									return (a)
									}
			)

allfam=do.call(rbind, fams)


###########################################
##   Assign families and superfamilies   ##
###########################################

mcols(ltrs.mask)$name=paste(seqnames(ltrs.mask), start(ltrs.mask), end(ltrs.mask), sep='_')

## add the family column - different for the original intact level and the subtracted ones
mcols(ltrs.mask)$family=NA
mcols(ltrs.mask)$family[mcols(ltrs.mask)$nestlevel==0]=mapvalues(mcols(ltrs.mask)$name[mcols(ltrs.mask)$nestlevel==0], from=a$V2, to=a$V1)
mcols(ltrs.mask)$family[mcols(ltrs.mask)$name %in% allfam$V1]=mapvalues(mcols(ltrs.mask)$name[mcols(ltrs.mask)$name %in% allfam$V1], from=allfam$V1, to=allfam$V2)


## assign mtec families
pm=read.table(paste('families/', GENOMENAME, '.TEDB.8080.searchglobal.toponly.out', sep=''), header=F)
pm$wicker3=str_split_fixed(pm$V1, '_', 3)[,1]
pm=pm[pm$wicker3 %in% c('RLC', 'RLG', 'RLX'),] ### make sure it's from the order ### might be nice to compare orders to mine called by protein order?
pm$family=str_split_fixed(pm$V1, '_', 3)[,2]   ### get the family name
### put family names in instead of other copies
pm.gd=merge(pm, a, by.x='V2', by.y='V2', all.x=T)


### assign an mtec family if any of the family members of mine matched
mcols(ltrs.mask)$mtec.fam=NA
for (i in 1:nrow(pm.gd)){
  mcols(ltrs.mask)$mtec.fam[mcols(ltrs.mask)$family==pm.gd$V1.y[i]]=pm.gd$family[i]
}

### for jeff
mcols(ltrs.mask)$famsize=table(mcols(ltrs.mask)$family)[mcols(ltrs.mask)$family]
mcols(ltrs.mask)$mtec.fam[is.na(mcols(ltrs.mask)$mtec.fam) & mcols(ltrs.mask)$famsize==1 & !is.na(mcols(ltrs.mask)$famsize)][2013]='alicia'


## for each copy, assign whether it can be attributed to either superfamily based on protein order
copia=c('GAG.AP.INT.RT.RNaseH.ENV', 'GAG.AP.INT.RT.RNaseH', 'INT.RT.RNaseH', 'AP.INT.RT.RNaseH', 'GAG.GAGCOAT.GAG.AP.INT.RT.RNaseH.ENV', 'GAG.AP.RT.RNaseH', 'AP.INT.RT.RNaseH.ENV', 'GAG.GAGCOAT.GAG.AP.INT.RT.RNaseH', 'INT.RT
.RNaseH.ENV')
gypsy=c('GAG.AP.RT.RNaseH.INT', 'GAG.AP.RT.RNaseH.INT.CHR', 'GAG.GAGCOAT.AP.RT.RNaseH.INT.CHR', 'RT.RNaseH.INT.CHR', 'GAG.GAGCOAT.GAG.AP.RT.RNaseH.INT.CHR', 'RT.RNaseH.INT', 'GAG.AP.RT.INT', 'AP.RT.RNaseH.INT.CHR', 'GAG.GAGCOAT.AP.RT.RNaseH.INT', 'GAG.RT.RNaseH.INT.CHR', 'RNaseH.INT', 'GAG.RT.RNaseH.INT', 'GAG.AP.RT.INT.CHR', 'RNaseH.INT.CHR', 'AP.RT.RNaseH.INT')
mcols(ltrs.mask)$superfam=NA
mcols(ltrs.mask)$superfam[prot.mask$geneorder %in% copia]='RLC'
mcols(ltrs.mask)$superfam[prot.mask$geneorder %in% gypsy]='RLG'


## assign superfamily to each family based on whether protein order can be assigned
fam.sup=sapply(names(table(mcols(ltrs.mask)$family)), function(x) names(which.max(table(mcols(ltrs.mask)$superfam[which(mcols(ltrs.mask)$family==x)]))))
mcols(ltrs.mask)$fam.sup=mapvalues(mcols(ltrs.mask)$family, from=names(unlist(fam.sup)), to=unlist(fam.sup))
mcols(ltrs.mask)$fam.sup[!mcols(ltrs.mask)$fam.sup %in% c('RLC', 'RLG') & !is.na(mcols(ltrs.mask)$famname)]='RLX'  ## only do for those that already have a family!!

### switch the LTR naming of families from silix to one that reflects superfamily ala Wicker et al. (2007), and is 5 digits long
mcols(ltrs.mask)$famname=sapply(1:length(ltrs.mask), function(x) sub('LTR0', mcols(ltrs.mask)$fam.sup[x], mcols(ltrs.mask)$family[x]))

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
for (x in names(table(mcols(ltrs.mask)$famname))){
  mcols(ltrs.mask)$Name[mcols(ltrs.mask)$famname==x & !is.na(mcols(ltrs.mask)$famname)]=paste(x, 'B73v4', str_pad(1:sum(mcols(ltrs.mask)$famname[!is.na(mcols(ltrs.mask)$famname)]==x), 5, pad='0'), sep='')
}


### this is the original name with positions on the subtracted chromosomes (not the coordinates from the full v4 chromosomes)
prot.mask$origname=paste(seqnames(ltrs.mask), prot.mask$element.start, prot.mask$element.end, sep='_')
mcols(ltrs.mask)$origname=prot.mask$origname


#### output a table so that I can recluster these. will need names for fasta headers and subtracted level
### write them separately for unplaced and chromosomes so that I can split them easily.
write.table(data.frame(prot.mask$nestlevel[is.na(mcols(ltrs.mask)$family) & prot.mask$sequence %in% as.character(1:10)], prot.mask$origname[is.na(mcols(ltrs.mask)$family)& prot.mask$sequence %in% as.character(1:10)]), 'copies_without_family_chromosome.txt', quote=F, row.names=F, col.names=F, sep='\t')
write.table(data.frame(prot.mask$nestlevel[is.na(mcols(ltrs.mask)$family) & !prot.mask$sequence %in% as.character(1:10)], prot.mask$origname[is.na(mcols(ltrs.mask)$family)& !prot.mask$sequence %in% as.character(1:10)]), 'copies_without_family_unplaced.txt', quote=F, row.names=F, col.names=F, sep='\t')
## this is for regular, all in one
##write.table(data.frame(prot.mask$nestlevel[is.na(mcols(ltrs.mask)$family)], prot.mask$origname[is.na(mcols(ltrs.mask)$family)]), 'copies_without_family.txt', quote=F, row.names=F, col.names=F, sep='\t')


### save these to update once i have families for all!
saveRDS(prot.mask, paste(GENOMENAME, '_prot.mask.dataframe.RDS', sep=''))
saveRDS(ltrs.mask, paste(GENOMENAME, '_ltrs.mask.genomicranges.RDS', sep=''))




#### gff trial
ltrs.gff=data.frame(seqnames(ltrs.mask), 'LTRharvest', 'LTR_retrotransposon', start(ltrs.mask), end(ltrs.mask), '.', strand(ltrs.mask), '.', paste('ID=', mcols(ltrs.mask)$Name, ';Name=', paste(mcols(ltrs.mask)$Name, mcols(ltrs.mask)$mtec.fam, paste('LTRsimilarity', mcols(ltrs.mask)$ltr_similarity, sep=''), sep='_'), sep=''))
write.table(ltrs.gff, paste(GENOMENAME, '.nonsubtractedLTR.lackingfamilies.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)




###### need to run cluster_distrupted_without-families.sh in the families directory, to get families assigned for each of these!

nf=read.table(paste('families/', GENOMENAME, '.nofamily.5ltr.8080.fnodes', sep=''))
nf$V1=as.character(nf$V1)
nf$V2=as.character(nf$V2)

nofam=is.na(mcols(ltrs.mask)$family)


mcols(ltrs.mask)$family[nofam]=mapvalues(mcols(ltrs.mask)$origname[nofam], from=nf$V2, to=nf$V1)
mcols(ltrs.mask)$superfam[nofam & prot.mask$geneorder %in% copia]='RLC'
mcols(ltrs.mask)$superfam[nofam & prot.mask$geneorder %in% gypsy]='RLG'


## assign superfamily to each family based on whether protein order can be assigned
fam.sup=sapply(names(table(mcols(ltrs.mask)$family)), function(x) names(which.max(table(mcols(ltrs.mask)$superfam[which(mcols(ltrs.mask)$family==x)]))))
mcols(ltrs.mask)$fam.sup[nofam]=mapvalues(mcols(ltrs.mask)$family[nofam], from=names(unlist(fam.sup)), to=unlist(fam.sup))
mcols(ltrs.mask)$fam.sup[!mcols(ltrs.mask)$fam.sup %in% c('RLC', 'RLG') & nofam]='RLX'  ## only do for those that already have a family!!


#########NEEEEED TO RENUMBER THESE!!!!  ################# right now just putting a 1 in front
### switch the LTR naming of families from silix to one that reflects superfamily ala Wicker et al. (2007), and is 5 digits long
mcols(ltrs.mask)$famname[nofam]=sapply(1:sum(nofam), function(x) sub('DLTR00', paste(mcols(ltrs.mask)$fam.sup[nofam][x], 1, sep=''), mcols(ltrs.mask)$family[nofam][x]))

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
for (x in names(table(mcols(ltrs.mask)$famname[nofam]))){
  mcols(ltrs.mask)$Name[mcols(ltrs.mask)$famname==x & nofam]=paste(x, 'B73v4', str_pad(1:sum(mcols(ltrs.mask)$famname[nofam]==x), 5, pad='0'), sep='')
}


## final gff

ltrs.gff=data.frame(seqnames(ltrs.mask), 'LTRharvest', 'LTR_retrotransposon', start(ltrs.mask), end(ltrs.mask), '.', strand(ltrs.mask), '.', paste('ID=', mcols(ltrs.mask)$Name, ';Name=', paste(mcols(ltrs.mask)$Name, mcols(ltrs.mask)$mtec.fam, paste('LTRsimilarity', mcols(ltrs.mask)$ltr_similarity, sep=''), sep='_'), sep=''))
write.table(ltrs.gff, paste(GENOMENAME, '.allFamilies.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

saveRDS(ltrs.mask, paste(GENOMENAME, '.ltrs.mask.GRanges.RDS', sep=''))

#############################################
###    GFF WITH GAPS FOR DISRUPTED TES    ###
#############################################


### generate a gff with 

ltrs.dj=disjoin(ltrs.mask, ignore.strand=T)
ltrs.overlap=findOverlaps(ltrs.dj, ltrs.mask)
mcols(ltrs.dj)=splitAsList(mcols(ltrs.mask)$name[subjectHits(ltrs.overlap)],queryHits(ltrs.overlap))
#mcols(ltrs.dj)$ltrsim=splitAsList(mcols(ltrs.mask)$ltr_similarity[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$ltrsim=splitAsList(mcols(ltrs.mask)$ltr_similarity[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$ltrsim.c=sapply(mcols(ltrs.dj)$ltrsim, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$ltrsim.self=sapply(mcols(ltrs.dj)$ltrsim.c, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])

mcols(ltrs.dj)$geneorder=splitAsList(prot.mask$geneorder[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$geneorder.c=sapply(mcols(ltrs.dj)$geneorder, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$geneorder.self=sapply(mcols(ltrs.dj)$geneorder.c, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])


mcols(ltrs.dj)$nestlevel=splitAsList(mcols(ltrs.mask)$nestlevel[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$nestlevel=sapply(mcols(ltrs.dj)$nestlevel, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$nestlevel=sapply(mcols(ltrs.dj)$nestlevel, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])



mcols(ltrs.dj)$ownname=splitAsList(mcols(ltrs.mask)$Name[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$ownname=sapply(mcols(ltrs.dj)$ownname, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$ownname=sapply(mcols(ltrs.dj)$ownname, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])


mcols(ltrs.dj)$famname=splitAsList(mcols(ltrs.mask)$famname[subjectHits(ltrs.overlap)], queryHits(ltrs.overlap))
mcols(ltrs.dj)$famname=sapply(mcols(ltrs.dj)$famname, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$famname=sapply(mcols(ltrs.dj)$famname, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])


ltrs.overlap.equal=findOverlaps(ltrs.dj, ltrs.mask, type='equal')
ltrs.gff=data.frame(seqnames(ltrs.dj), '.', 'region', start(ltrs.dj), end(ltrs.dj), '.', strand(ltrs.dj), '.', 'col9')
names(ltrs.gff)=c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')
ltrs.gff$V3=as.character(ltrs.gff$V3)
ltrs.gff$V9=as.character(ltrs.gff$V9)
mcols(ltrs.dj)$value=sapply(mcols(ltrs.dj)$value, function(x) paste(as.character(x), collapse=','))
mcols(ltrs.dj)$nestlevels=str_count(mcols(ltrs.dj)$value, ',')+1
for (linenum in 1:length(ltrs.dj)){
  if ( linenum %in% queryHits(ltrs.overlap.equal)){
    ltrs.gff[linenum,3]='LTR_retrotransposon'
    ltrs.gff[linenum,9]=paste('ID=', mcols(ltrs.dj)$ownname[linenum], 
#    					'; Parent=', mcols(ltrs.dj)$value[linenum], 
    					'; Name=', mcols(ltrs.dj)$famname[linenum],'_LTRsimilarity', mcols(ltrs.dj)$ltrsim.self[linenum], '_GENEORDER-', mcols(ltrs.dj)$geneorder.self[linenum], ',overlappingTEcopies', mcols(ltrs.dj)$nestlevel[linenum]
    					, sep='')
#    					,
						
  }
  else{
    ltrs.gff[linenum,3]='transposon_fragment'
    ltrs.gff[linenum,9]=paste('ID=FRAG_', mcols(ltrs.dj)$ownname[linenum], 
#    					'; Parent=', mcols(ltrs.dj)$value[linenum], 
    					'; Name=,', mcols(ltrs.dj)$famname[linenum],'_LTRsimilarity', mcols(ltrs.dj)$ltrsim.self[linenum], '_GENEORDER-', mcols(ltrs.dj)$geneorder.self[linenum], ',overlappingTEcopies', mcols(ltrs.dj)$nestlevel[linenum]
    					, sep='')
    					
  }
}




write.table(ltrs.gff, paste(GENOMENAME, '.nested_subtracted_ltr.family.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)


