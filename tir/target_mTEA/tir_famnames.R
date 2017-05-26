

tir=tir[-rmRows,]
tir.gr=tir.gr[-rmRows,]

tir.gr$mtec=tir$mtec
tir.gr$mtecfamnum=substr(tir.gr$mtec, 7,11)
tir.gr$famname=paste(tir.gr$superfam, tir.gr$mtecfamnum, sep='')

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
mcols(tir.gr)$Name=NULL
for (x in names(table(mcols(tir.gr)$famname))){
  mcols(tir.gr)$Name[mcols(tir.gr)$famname==x & !is.na(mcols(tir.gr)$famname)]=paste('B73v4', str_pad(1:sum(mcols(tir.gr)$famname[!is.na(mcols(tir.gr)$famname)]==x), 5, pad='0'),sep='')
}

tir.gr$sup=substr(tir.gr$mtec, 1,3)
tir.gr$ID=paste0(tir.gr$sup, tir.gr$Name)


te=readRDS('../../highconf_filtered_te_set_final.RDS')
newtir=findOverlaps(te, tir.gr, type='equal')
te$IDnew=te$ID
te$IDnew[queryHits(newtir)]=tir.gr$ID[subjectHits(newtir)]

te$Namenew=te$Name
te$Namenew[queryHits(newtir)]=paste(tir.gr$ID[subjectHits(newtir)], paste0('TSDlen', tir$tsdlen)[subjectHits(newtir)], paste0('TIRlen', tir$tirlength)[subjectHits(newtir)], sep='_')

te$ID=te$IDnew
te$IDnew=NULL

te$Name=te$Namenew
te$Namenew=NULL

names(table(te$ID))[table(te$ID)==2]

export.gff3(te, 'B73v4_structrural_filtered_newTIRID.Oct2516.gff3')



#NOW, remove the SINE dups!

sinefree=te[-queryHits(findOverlaps(te, type='equal', drop.self=T))[c(T,F)],]
export.gff3(sinefree, 'B73v4_structrural_filtered_newTIRID_noSINEdup.Oct2516.gff3')
