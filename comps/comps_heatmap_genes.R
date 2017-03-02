library("ggplot2");
library("reshape2");
library("RColorBrewer");
library("amap")

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args
input_fname = args[1];
output_fname = args[2];
num_control = as.numeric(args[3]);
num_test = as.numeric(args[4]);
comps_id = args[5];
offset = 10

df <- read.csv(file=input_fname, header=T, quote="$", sep="\t");

#Subset
subset.df <- df[1:100, ];

subset.data <- as.matrix(subset.df[, 10:(10+num_control+num_test-1)]);
#norm.data <- log2(subset.data / rowMeans(subset.data));
norm.data <- subset.data / rowMeans(subset.data);
rownames(norm.data) <- subset.df$gene_name;

norm.dist <- Dist(norm.data, method="pearson");

norm.df <- subset.df;
norm.df[, 10:(10+num_control+num_test-1)] <- norm.data;

#Cluster
clust.df <- hclust(norm.dist, method="single");

#Order
order.df <- norm.df[clust.df$order, ];
order.df$gene_name <- factor(order.df$gene_name, levels=order.df$gene_name[clust.df$order]);

melt.df <- melt(order.df, measure.vars=colnames(df)[10:(10+num_control+num_test-1)]);
melt.df$gene_name <- factor(melt.df$gene_name, levels=order.df$gene_name)

#Plot
pdf(file=paste(output_fname, ".pdf", sep=""), width=6, height=18);
# scale_colour_gradient(trans="log")
# scale_fill_gradientn(colours=brewer.pal(11,"RdBu"))
#+ ggtitle("Heatmap")
curr.plot <- ggplot(melt.df, aes(x=variable, y=gene_name, fill=value)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_tile() + scale_fill_gradientn(colours=brewer.pal(3,"RdBu")) + xlab("") + ylab("") + theme(plot.title=element_text(size=10)) + ggtitle(comps_id)
print(curr.plot);
#plot(clust.df)
dev.off();

png(file=paste(output_fname, ".png", sep=""), width=600, height=1500);
curr.plot <- ggplot(melt.df, aes(x=variable, y=gene_name, fill=value)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_tile() + scale_fill_gradientn(colours=brewer.pal(3,"RdBu")) + xlab("") + ylab("") + theme(plot.title=element_text(size=10)) + ggtitle(comps_id)
print(curr.plot);
dev.off();
