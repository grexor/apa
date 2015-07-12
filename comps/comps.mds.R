library(edgeR)

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args
input_filename = args[1];
output_filename = args[2];
num_control = as.numeric(args[3]);
num_test = as.numeric(args[4]);
offset = 10

data <- read.delim(input_filename,row.names=NULL)
data = data[, offset:(offset+num_control+num_test-1)]

require(graphics)
svg(output_filename, width=10, height=7, pointsize=9, bg="transparent")
op <- par(mar = par("mar") + c(1,1,1,8))
plotMDS(data)
dev.off()
