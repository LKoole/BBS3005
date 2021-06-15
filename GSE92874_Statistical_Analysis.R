# Clear Global Environment
rm(list=ls())

# Set working directory
setwd("/Users/lisakoole/Desktop/Datasets/GSE92874")

#import data from file
file <- file.path(getwd(),"GSE92874_RNA.txt")
dat <- read.delim(file, header = TRUE, sep = "\t", dec = ".")


## Format the data ##
# Create experimental group factor variable
desc <- cbind(group=c(rep("control",4),rep("schiz", 4)))
desc <- as.data.frame(cbind(desc,individual=paste("s", rep(1:4, 2), sep = "")))
rownames(desc) <- paste(desc$group,desc$individual,sep="_")

# Add sample and gene names to data.frame
colnames(dat) <- paste(c("gene", rownames(desc), "log2_fold_change", 
                         "p_value", "q_value", "significant"))
rownames(dat) <- paste("gene",1:dim(dat)[1],sep="")

# Load needed libraries
library(org.Hs.eg.db)
library(dplyr)

# In case of duplicates, remove the least significant ones
dat.ord <- dat[order(dat$p_value),]
dat2 <- dat.ord %>% distinct(gene, .keep_all = TRUE)

## Add annotation ##
columns(org.Hs.eg.db)
ann <- AnnotationDbi::select(org.Hs.eg.db,keys=dat.ord[,1], keytype = "SYMBOL", 
                             columns=c("ENTREZID", "SYMBOL"))

# Make sure there are no duplicates regarding the identifier 
# that will be used as identifier
ann.noRep <- ann %>% distinct(SYMBOL, .keep_all = TRUE)

# Add annotations to data.frame
table(ann.noRep$SYMBOL==dat2[,1])
dat2 <- cbind(dat2, entrezid=ann.noRep[,2])

# Export data to txt files 
write.table(dat2, "GSE92874_annotated_noRep.txt", sep = "\t", dec = ".", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
write.table(dat2[,c(1, 10, 11, 12, 14)], "GSE92874_for_pathvisio_.txt", sep = "\t", quote = FALSE,
            dec = ".", row.names = FALSE, col.names = TRUE)


## Create factors and countdata ##
factors <- c("group", "individual")
countdata <- dat2[,c(2:9)]
rownames(countdata) <- paste(dat2[,1])


## Filter lowely expressed genes ##
thresh <- countdata > 0

# Summary of how many TRUEs there are in each row. 
table(rowSums(thresh))
# In this case:  15268 are expressed in all 8 samples, 3238 are not expressed in any sample

# Keep genes that are expressed in all samples
keep <- rowSums(thresh) == 8 
summary(keep)

countdata.expr <- countdata[keep,]
dat.expr <- dat2[keep, ]

write.table(dat.expr, "GSE92874_filtered_data.txt", sep = "\t", dec = ".", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(dat.expr[,c(1,10,11,12,14,15)], "GSE92874_filtered_data2.txt", sep = "\t", dec = ".", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Transform counts to log2 scale 
logcounts <- log(countdata.expr, 2)
logcounts <- do.call(data.frame,                      # Replace Inf in data by NA
                     lapply(logcounts,
                            function(x) replace(x, is.infinite(x), NA)))
rownames(logcounts) <- paste(rownames(countdata.expr))


## Create plots ##
# Load needed libraries
library(Glimma)
library(RColorBrewer)
library(NMF)
library(gplots)
library(ggplot2)
library(ggforce)

# Create boxplot of logcounts
boxplot(logcounts, xlab = "", ylab="log2(FPKM)", las=2)
# Calculate medians per sample
logcounts.median <- as.data.frame(lapply(logcounts, median, na.rm = TRUE))
# Calculate median for all data
logcounts.median.all <- median(as.matrix(logcounts))

# Add blue lines corresponding to the median of logcounts
abline(h=logcounts.median.all,col="blue")
title("Boxplots of logFPKMs (outliers not removed)")


## Determine outliers ##
outliers <- boxplot(logcounts, plot=FALSE)$out
logcounts.noOut <- logcounts[-which(as.matrix(logcounts) %in% outliers),]

# Create Boxplot with outliers 
boxplot(logcounts.noOut, xlab="", ylab="log2(FPKM)", las=2)
title("Boxplots of logFPKMs (outliers removed)")
# In this case: not included in the results 
dev.off()

## Create MDS plot ##
labels <- paste(desc$group, desc$individual)
glMDSPlot(logcounts, labels=labels, groups=factor(desc$group), folder="mds")

## Create separate PCA plot ## 
library(ggfortify)
pca_res <- prcomp(t(logcounts), scale. = TRUE)
# x is the PCA number on x axis, y is PCA number on y-axis
autoplot(pca_res, labels=labels, data = t(logcounts), colour = col.cell,
         label = FALSE, size = 7, x = 2, y = 3)

## Hierarchical clustering with heatmap ##
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset countdata matrix
highly_variable_lfpkm <- as.matrix(logcounts[select_var,])
# Alternatively, use q values to determine most variables ones
## highly_variable_lfpkm <- as.matrix(dat.expr[1:500,c(2:9)])
## rownames(highly_variable_lfpkm) <- dat.expr[1:500,1]
dim(highly_variable_lfpkm)
head(highly_variable_lfpkm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("dimgrey","grey")[as.factor(desc$group)]

# Plot the heatmap
heatmap.2(highly_variable_lfpkm, col=morecols(50),trace="none", 
          main="Top 500 most variable genes across samples",
          ColSideColors=col.cell, scale="row",
          density.info = "none",
          margins =c(7,7))
dev.off()


# CLUSTERING #
# Load needed libraries
library(bioDist)
library(factoextra)

# Correlation-based distance method
res.dist <- get_dist(t(logcounts), method = "pearson")
hc = hclust(res.dist, method = "ward.D2")
plot(hc, main = "Clustering of brain tissue genes")

# Visualize the dissimilarity matrix
# Note: Red color corresponds to small distance and 
# blue color indicates big distance between observation
fviz_dist(res.dist, lab_size = 8)
# Heatmap of sample groups with cluster dendrograms
heatmap.2(as.matrix(res.dist), col = morecols(50), scale = "none", symm = TRUE, 
          trace = "none", density.info = "none", dendrogram = "column",
          margins=c(7,7), RowSideColors=col.cell)
dev.off()


## Determine differentially expressed genes ##
diff.genes <- filter(dat.expr, q_value < 0.05 & (log2_fold_change > 1 | log2_fold_change < -1))
# Eliminate NA values
diff.genes.noNA <- diff.genes[which(!diff.genes$entrezid=="NA"),]

# Export list of differentially expressed genes to txt files
write.table(diff.genes.noNA[,1], "GSE92874_DEG_filtered.txt", sep = "\t", dec = ".", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(diff.genes.noNA[,15], "GSE92874_DEG_filtered_UP.txt", sep = "\t", dec = ".", 
            quote = FALSE, row.names = FALSE, col.names = "UniprotID")
write.table(dat2[,1], "GSE92874_background_genes.txt", sep = "\t", dec = ".", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

                            
## GO annotation of DEGs ##
annGO.Deg <- AnnotationDbi::select(org.Hs.eg.db,keys=diff.genes[,1], keytype = "SYMBOL", 
                                columns=c("ENTREZID", "SYMBOL", "GO"))
# Extract named vector of all terms
library(GO.db)
goterms <- Term(annGO.Deg[,3])
annGO.Deg <- as.data.frame(cbind(annGO.Deg,Terms=paste(goterms)))
# Add terms to data frame
annGO.Deg.BP <- annGO.Deg[which(annGO.Deg$ONTOLOGY=="BP"),]

write.table(annGO.Deg.BP, "GSE92874_DEG_GOBP.txt", sep = "\t", dec = ".", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
    
                            
## Create MD plots ##
par(mfrow = c(1,2))
plotMD(logcounts, column = 5)
abline(h=0,col="grey")
plotMD(logcounts, column = 6)
abline(h=0,col="grey")


## Volcano plot ##
deseq.threshold <- as.factor((dat.expr$log2_fold_change < -1 | dat.expr$log2_fold_change > 1) & 
  dat.expr$q_value < 0.05)
deseq.upreg <- as.factor(dat.expr$log2_fold_change > 1 & dat.expr$q_value < 0.05)
deseq.downreg <- as.factor(dat.expr$log2_fold_change < -1 & dat.expr$q_value < 0.05)
deseq <- as.factor(paste(deseq.downreg, deseq.upreg))

# Get the indices of interesting genes using the threshold
xi <- which(deseq.threshold == TRUE)

# Number of interesting genes
length(xi)

# Load ggrepel package
library(ggrepel)

# Set up colour vector for celltype variable
col.deg <- c("gray25", "indianred", "cornflowerblue")[deseq]
fc <- log2_fold_change
padj <- q_value

g = ggplot(data=dat.expr, 
           aes(x=log2_fold_change, y=-log10(q_value))) +
  # Define the look of the points
  geom_point(alpha=0.4, size=1.75, colour=col.deg) +
  # Hide the legend
  theme(legend.position = "none") +
  # Apply another theme
  theme_bw() + theme(legend.position="none") +
  # Add the lines separating the DEGs
  geom_vline(xintercept = 1, color = 'gray69') +
  geom_vline(xintercept = -1, color = 'gray69') +
  geom_hline(yintercept = -log10(0.05), color = 'gray69') +
  xlab("Log2 fold change (Log2FC)") + ylab("-log10(Padj)") +
  ggtitle("IPSC trisomic vs disomic DEGs according to cuffdiff\nwith FDR <= 0.05 and absolute FC >= 1 or FC <= -1") +
  # Add the text labels only for DEGs. Using `ggrepel` to make sure labels don't overlap
  # Note: p-value threshold is set very low, to only highlight very significant genes
  geom_text_repel(aes(label=ifelse(q_value < 0.1e-2 & (log2_fold_change < -1 | log2_fold_change > 1),
                                   gene, '')))

# Show the plot
g
dev.off()

## Alternative method of making a volcano plot ##
# Base plot
with(dat.expr, plot(fc, -log10(padj), pch=20,
                         main=paste("IPSC trisomic vs disomic DEGs according to cuffdiff",
                                    "with FDR <= 0.05 and absolute FC >= 1", sep="\n"), xaxt = 'n'))

# Add custom x-axis with more ticks (default was (-)2, 4, 6)
# NOTE: adjust for your own data!
axis(side = 1, at = seq(-6, 9, 1))

with(subset(dat.expr, q_value < 0.05), 
     # Add points for statistically significant genes
     points(log2_fold_change, -log10(q_value), pch=20, col="red"))

with(subset(dat.expr, log2_fold_change > 1 | log2_fold_change < -1), 
     # Add points for expression > logfc_threshold
     points(log2_fold_change, -log10(q_value), pch=20, col="orange"))

with(subset(dat.expr, q_value < 0.05 & abs(log2_fold_change) > 1), 
     # Add points based on both significance and expression
     points(log2_fold_change, -log10(q_value), pch=20, col="green"))

# Add threshold lines
abline(h=-log10(0.05), lwd=1.5, lty=2, col='gray')
abline(v=c(-1, 1), lwd=1.5, lty=2, col='gray')
dev.off()
