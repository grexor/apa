library("edgeR")
library("data.table")

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args
input_fname = args[1];
output_fname = args[2];
num_control = as.numeric(args[3]);
num_test = as.numeric(args[4]);
offset = 10

gx <- fread(input_fname, colClasses=list(character=1:1))
gxcounts <- gx[,-c(1:9),with=FALSE]

group <- factor(c(rep(c(1), num_control), rep(c(2), num_test)))
d <- DGEList(counts=gxcounts, group=group)
d <- calcNormFactors(d)

d$genes <- gx[,1:9,with=FALSE]
d$genes <- cbind(d$genes, round(cpm(d),1))

d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
d <- exactTest(d)
d = topTags(d, n=Inf)

write.table(d, file=output_fname, sep="\t", row.names=FALSE, quote=FALSE)
