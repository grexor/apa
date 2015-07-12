library("edgeR")
library("data.table")

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args
fname = args[1];
fbase = sub("^([^.]*).*", "\\1", fname) 

#------------
# gene level
#------------

gx <- fread(fname, colClasses=list(character=1:1))
gxcounts <- gx[,-c(1:9),with=FALSE]

d <- DGEList(gxcounts)
d <- calcNormFactors(d)
d$genes <- gx[,1:9,with=FALSE]
cps <- cpm(d)
d$genes <- cbind(d$genes, round(cps,1))

write.table(d$genes, file=paste(fname, ".genes_edgeR.tab", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=colnames(d$genes))

svg(paste(fbase, ".genes_mds.svg", sep=""), width=10, height=7, pointsize=9, bg="transparent")
plotMDS(d)
dev.off()


#------------
# site level
#------------

#fname = args[2];
#fbase = sub("^([^.]*).*", "\\1", fname) 
#gx <- fread(paste(fbase, ".sites.tab", sep=""), colClasses=list(character=1:1))
#gxcounts <- gx[,-c(1:6),with=FALSE]
#d <- DGEList(gxcounts)
#d <- calcNormFactors(d)
#d$genes <- gx[,1:6,with=FALSE]
#cps <- cpm(d)
#d$genes <- cbind(d$genes, round(cps,1))
#write.table(d$genes, file=paste(fbase, ".sites_cpm.tab", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=colnames(d$genes))
