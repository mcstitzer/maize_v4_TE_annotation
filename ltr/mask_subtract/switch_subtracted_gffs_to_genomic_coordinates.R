library(rtracklayer)
## want to use lots of thresholds
maxSubtract=80
maxSubtract=as.numeric(commandArgs(trailingOnly=TRUE))

GENOMENAME='B73V4.pseudomolecule'

### use the samtools faidx generated index to know the length of each contig/chromosome
gs=read.table(paste('../../', GENOMENAME, '.fasta.fai', sep=''), sep='\t', header=F)

### gl is a list of vectors (initalized with NA) that are the length of each contig/chromosome
###    later, I will update each vector with the position at which a TE exists and is subtracted
gl=lapply(gs$V2, function(x) rep(NA, x))
names(gl)=gs$V1     ### here, I'm naming the lists so they can be accessed by contig/chromosome name
### index is a list of vectors, will be updated to keep track of unvisited bases (those that are in the subtracted chromosome/contig fasta)
index=lapply(gl, function(x) is.na(x))  

############## 
## get all the file names ready
##############
##!!!!!!!!!! MAKE SURE THESE ARE THE LTRRETROTRANSPOSON RECORD ONLY GFFS!!!!!! WAS HAVING LOTS OF TROUBLE WITH IMPROPER COORDINATES UNTIL I FIXED THIS

masked=paste('subtract', 1:maxSubtract, '/', GENOMENAME, '.subtract', 1:maxSubtract, '.ltrharvest.contignames.tsd.ltrretrotransposon.gff3', sep='')
filenames=c(paste(GENOMENAME, '.ltrharvest.contignames.tsd.ltrretrotransposon.gff3', sep=''), masked)


##############
## loop through each gff, by chromosome/contig, updating positions
##############

ranges=c()  ### initialize this value, it will represent the bp coordinates we need to mask
for (x in 1:(maxSubtract+1)){stageNum=x
                gff=import.gff(filenames[x], version='3', asRangedData=F)  ### import a gff
                index=lapply(gl, function(j) is.na(j))       ### update index (which bases are unvisited (still NA))
                for (i in seqlevels(gff)){                   #### go by chromsome/contig
					g=gff[seqnames(gff)==i]     ### get only those entries from this chromosome/contig
					ranges=unlist(lapply(g, function(r) start(r):(end(r)))) ### record every position that is present in this annotation, in the coordinates of the subtracted fasta
					gl[[i]][which(index[[i]])[ranges]]=stageNum             ### update the positions in gl with their subtraction level
					}
				ranges=c()   ### clear this variable
				rm(gff)      ### was having memory troubles, so remove the imported file, and garbage collect
				gc()
				}                
saveRDS(gl, paste(GENOMENAME, '_gl_', maxSubtract, '.RDS', sep=''))

#################################
###   REWRITE GFFS          #####
#################################
# replacing the subtracted coordinates with the real coordinates from the chromosomes/contigs

### use different file names because I want these to include stuff like ltr position and tsd position
###   could also use to switch ltrdigest results, if desired.

gffmasked=paste('subtract', 1:maxSubtract, '/', GENOMENAME, '.subtract', 1:maxSubtract, '.ltrharvest.contignames.gff3', sep='')
gfffilenames=c(paste('../', GENOMENAME, '.ltrharvest.contignames.gff3', sep=''), gffmasked)

for (sublev in 1:length(gfffilenames)){
	gff=import.gff(gfffilenames[sublev], version='3', asRangedData=F)
	for (i in seqlevels(gff)){
		rowindex=as.vector(seqnames(gff)==i)   ## only want to change positions in the gff for this chromosome/contig
		chrlist=gl[[i]]    ### these are the gff number at which each bp was annotated as TE
		chrlist.subtract=which(is.na(chrlist) | chrlist>=sublev)   ### these are the positions that would have been visible to ltrharvest when annotating this TE
		end(gff)[rowindex]=chrlist.subtract[end(gff)[rowindex]]   ### do end first so genomic ranges doesn't complain about negative widths
		start(gff)[rowindex]=chrlist.subtract[start(gff)[rowindex]]
		}
	export.gff(gff, paste(gfffilenames[sublev], '.contigpositions.gff3', sep=''), 'gff3')
	}


