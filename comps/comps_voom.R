options(digits=3, width=120)
library(limma)
library(edgeR)

args = commandArgs(trailingOnly = T);
folder = args[1];
comps_id = args[2];
input_filename = args[3];
num_control = as.numeric(args[4]);
num_test = as.numeric(args[5]);
offset = 9

data <- read.delim(input_filename)
genes = data[,4]
sites = data[,7]

data = data[,offset:(offset+num_control+num_test-1)]
colnames(data)

group <- factor(c(rep(c(1), num_control), rep(c(2), num_test)))
group

d = DGEList(counts=data, group=group)
d = calcNormFactors(d, method="TMM")

design = model.matrix(~factor(c(rep(c(1), num_control), rep(c(2), num_test))))
colnames(design) = c("intercept", "group")
design

bitmap(file=paste(folder, "/", comps_id, ".voom_mv.png", sep=""), width=1024, height=768, type="png16m", units="px", pointsize=16, taa = 4, gaa = 4)
y = voom(d, design, plot=T)
dev.off()

svg(paste(folder, "/", comps_id, ".voom_mds.svg", sep=""), bg = "transparent")
plotMDS(y)
dev.off()

fit = lmFit(y, design=design)
fit.de = eBayes(fit, robust=T)
summary(decideTests(fit.de))

ncol(fit)

colnames(fit)

ex = diffSplice(fit, geneid=genes, exonid=sites)
summary(decideTests(ex))

r.genes = topSplice(ex, coef=2, test="simes", n=Inf) # coef=2, compute change for test (1=control, 2=test); logFC >0: test > control
r.exons = topSplice(ex, coef=2, test="t", n=Inf)

write.table(r.genes, file=paste(folder, "/", comps_id, ".voom_genes.tab", sep=""), sep="\t", col.names=NA, quote=FALSE)
write.table(r.exons, file=paste(folder, "/", comps_id, ".voom_exons.tab", sep=""), sep="\t", col.names=NA, quote=FALSE)

# skip SVGs for 10 top genes
#for (i in 1:10)
#{
#    svg(paste(folder, "/", comps_id, ".i", i, ".svg", sep=""), bg = "transparent")
#    plotSplice(ex, rank = i)
#}

#svg(paste(folder, "/", comps_id, "ENSG00000067208.svg", sep=""), bg = "transparent")
#plotSplice(ex, geneid="ENSG00000067208")
#dev.off()
