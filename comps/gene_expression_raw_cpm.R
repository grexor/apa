library("edgeR")
library("data.table")

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args
input_fname = args[1]
output_fname = args[2]

gx <- fread(input_fname, colClasses=list(character=1:1))
gxcounts <- gx[,-c(1:2),with=FALSE]
d <- DGEList(counts=gxcounts)
d <- calcNormFactors(d)

d$genes <- gx[,1:2,with=FALSE]
d$genes <- cbind(d$genes, round(cpm(d), 1))

write.table(d$genes, file=output_fname, sep="\t", row.names=FALSE, quote=FALSE)
