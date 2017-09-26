library(data.table)
library(BiocParallel)

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args
input_folder = args[1];
output_fname = args[2];
num_control = as.numeric(args[3]);
num_test = as.numeric(args[4]);
comps_id = args[5];

BPPARAM = MulticoreParam(workers=12)
inDir = system.file(input_folder)

countFiles = c()
row_names = c()
conds = c()
libs = c()
for(i in 1:num_control) {
  countFiles = append(countFiles, paste(input_folder, "/", "dex_", comps_id, "_c", i, ".tab", sep=""))
  row_names = append(row_names, paste("c", i, sep=""))
  conds = append(conds, "control")
  libs = append(libs, "single-end")
}
for(i in 1:num_test) {
  countFiles = append(countFiles, paste(input_folder, "/", "dex_", comps_id, "_t", i, ".tab", sep=""))
  row_names = append(row_names, paste("t", i, sep=""))
  conds = append(conds, "test")
  libs = append(libs, "single-end")
}

sampleTable = data.frame(row.names = row_names, condition = conds, libType = libs)
suppressPackageStartupMessages( library( "DEXSeq" ) )
dxd = DEXSeqDataSetFromHTSeq(
countFiles,
sampleData=sampleTable,
design= ~ sample + exon + condition:exon)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )
dxd = testForDEU( dxd, BPPARAM=BPPARAM )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr1 = DEXSeqResults( dxd )
write.table(dxr1, file=output_fname, sep="\t", row.names=FALSE, quote=FALSE)
