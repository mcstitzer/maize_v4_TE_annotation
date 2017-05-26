library(rtracklayer)
## want to use lots of thresholds
maxSubtract=80
maxSubtract=as.numeric(commandArgs(trailingOnly=TRUE))

GENOMENAME='B73V4.pseudomolecule'

### use the samtools faidx generated index to know the length of each contig/chromosome
gs=read.table(paste('../../', GENOMENAME, '.fasta.fai', sep=''), sep='\t', header=F)

### gl is a list of vectors (initalized with NA) that are the length of each contig/chromosome
###    later, I will update each vector with the position at which a TE exists and is subtracted

gl=readRDS(paste(GENOMENAME, '_gl_', maxSubtract, '.RDS', sep=''))

#################################
###   REWRITE GFFS          #####
#################################
# replacing the subtracted coordinates with the real coordinates from the chromosomes/contigs

### use different file names because I want these to include stuff like ltr position and tsd position
###   could also use to switch ltrdigest results, if desired.

gffmasked=paste('subtract', 1:maxSubtract, '/', GENOMENAME, '.subtract', 1:maxSubtract, '.ltrharvest.contignames.gff3', sep='')
gfffilenames=c(paste('../', GENOMENAME, '.ltrharvest.contignames.gff3', sep=''), gffmasked)

for (sublev in 1:length(gfffilenames)){
	gff=import.gff(gfffilenames[sublev])
	for (i in seqlevels(gff)){
		rowindex=as.vector(seqnames(gff)==i)   ## only want to change positions in the gff for this chromosome/contig
		chrlist=gl[[i]]    ### these are the gff number at which each bp was annotated as TE
		chrlist.subtract=which(is.na(chrlist) | chrlist>=sublev)   ### these are the positions that would have been visible to ltrharvest when annotating this TE
		end(gff)[rowindex]=chrlist.subtract[end(gff)[rowindex]]   ### do end first so genomic ranges doesn't complain about negative widths
		start(gff)[rowindex]=chrlist.subtract[start(gff)[rowindex]]
		}
	export.gff(gff, paste(gfffilenames[sublev], '.contigpositions.gff3', sep=''), 'gff3')
	}


