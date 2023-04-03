# eLife-project (https-doi.org-10.7554-eLife.37892)
The R script is an exemplar script to explain spike-in normalisation in the context of ChIPseq 
#################### Import libraries  ####################################
library("devtools")
library("DESeq2")
library("Brundle")
library("BrundleData")
library("DESeq2")
library("ggplot2")
library("gplots")
library("dplyr")
library("ggrepel")
library("genefilter")
library("RColorBrewer")
library("GenomicAlignments")
library("BiocParallel")
library("vsn")
library("amap")
#################### Checking the quality of Chip-seq Count Data Data ####################################

#### Read the HTseq counts files 

setwd("/export/home1/users/mpb/amohamed/Mounted_dir/Project_mpb_pourlemomentnon/barneche/Ouardia/chip_seq_data/det1_project/chip_seq_diff_analysis/matrix_deseq2")

rawCountTable <- read.table("gene_count_matrix.bed")

cts <- as.matrix(rawCountTable)

#### Read the design file #############

coldata <- read.csv("design_file_det1_chip_seq_diff_analysis.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)
head(coldata)


############## DDS object ##################

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)

dds$group <- factor(paste0(dds$treatment, dds$genotype))
design(dds) <- ~ group

dds <- DESeq(dds, minReplicatesForReplace = Inf)

####Box plot of Data ####
par(mar=c(9,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, col="darkblue")

####Model distribution plot of Data ####
plotDispEsts(dds)

############## Regularized log transformation of Data ##################
rld <- rlog(dds)

############## Checking the quality of Data ##################

####Spearman correlation Hierarchical clustering ####
plot(hcluster(dist(t(counts(dds))), method= "spearman", link="average"), 
     cex=0.5, cex.main=1, main= "Hierarchical Clustering of det1-1 Chip-seq Count Data")

####PCA plot####
DESeq2::plotPCA(rld, intgroup = c("genotype", "treatment"))

####Heatmap plot based on Euclidean distance ####
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$genotype,
                                    rld$treatment,
                                    rld$replicate, sep="-")
colnames(sampleDistMatrix) <- paste(rld$genotype,
                                    rld$treatment,
                                    rld$replicate, sep="-")


colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(256)
heatmap.2(sampleDistMatrix, trace="none", col=colours, margins = c(9, 9)) 

####1000 most varying Gene expression Heatmap plot####
topVarGenest <- head(order(rowVars(assay(rld)), decreasing=TRUE), 1000)
heatmap.2(assay(rld)[topVarGenest, ], scale="row",labRow = FALSE,
          trace="none", dendrogram="column",
          col = colorRampPalette(c("green","black","red"))(256), margins = c(10, 10), 
          main="Sample Distance Matrix")


######################################################################
##########################Contrasts###################################
######################################################################



######################### Contrast WT_D_WT_L###############################################################


rawCountTable <- read.table("wt_d_wt_l_matrix.bed")

cts <- as.matrix(rawCountTable)


########################"DESeq2 matrix #####################################################################


coldata <- read.csv("design_file_wt_d_wt_l.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)
head(coldata)

############## DDS object ##################

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)

dds$treatment <- relevel(dds$treatment, ref = "D")
sizeFactors(dds) <- c(0.34, 0.35, 0.37, 0.32)
dds <- DESeq(dds, fitType='local')



######################### Contrast WT (D)  WT (L)###############################################################

res1 <- results(dds, alpha=0.05, cooksCutoff = FALSE, independentFiltering = FALSE)

table1 <- subset(res1, select=c("log2FoldChange", "padj", "pvalue"))
table1 <- as.data.frame(table1)
table1 <- na.omit(table1)
mcols(res1, use.names=TRUE)

ordered_res1 <- res1[order(res1$pvalue),]
summary(ordered_res1)

######################### Contrast det1-1 (D) det1-1 (L)####################################################


rawCountTable <- read.table("det_d_det_l_matrix.bed")

cts <- as.matrix(rawCountTable)


########################"DESeq2 matrix #####################################################################


coldata <- read.csv("design_file_det_d_det_l.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)
head(coldata)

############## DDS object ##################

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)

dds$treatment <- relevel(dds$treatment, ref = "L")
sizeFactors(dds) <- c(0.81, 0.37, 1.00, 0.91)
dds <- DESeq(dds, fitType='local')

res2 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)

table2 <- subset(res2, select=c("log2FoldChange", "padj", "pvalue"))
table2 <- as.data.frame(table2)
table2 <- na.omit(table2)
mcols(res2, use.names=TRUE)

ordered_res2 <- res2[order(res2$pvalue),]
summary(ordered_res2)

############## Contrast det1-1 (D) to WT (D) ###########

rawCountTable <- read.table("det_d_wt_d.bed")

cts <- as.matrix(rawCountTable)


########################"DESeq2 matrix #####################################################################


coldata <- read.csv("design_file_det_d_wt_d.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)
head(coldata)

############## DDS object ##################

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "WT")
sizeFactors(dds) <- c(0.81, 0.37,0.34, 0.35)
dds <- DESeq(dds, fitType='local')

res3 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)
table3 <- subset(res3, select=c("log2FoldChange", "padj", "pvalue"))
table3 <- as.data.frame(table3)
table3 <- na.omit(table3)
mcols(res3, use.names=TRUE)

ordered_res3 <- res3[order(res3$pvalue),]
summary(ordered_res3)


############## Contrast det1-1 (D) to WT (L) ###########


rawCountTable <- read.table("det_d_wt_l.bed")

cts <- as.matrix(rawCountTable)

########################"DESeq2 matrix #####################################################################


coldata <- read.csv("design_file_det_d_wt_l.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)

############## DDS object ##################

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)

dds$treatment <- relevel(dds$treatment, ref = "L")
sizeFactors(dds) <- c(0.81, 0.37,0.37, 0.32)
dds <- DESeq(dds, fitType='local')


res4 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE) 
table4 <- subset(res4, select=c("log2FoldChange", "padj", "pvalue"))
table4 <- as.data.frame(table4)
table4 <- na.omit(table4)
mcols(res4, use.names=TRUE)

ordered_res4 <- res4[order(res4$pvalue),]
summary(ordered_res4)


############## Contrast det1-1 (L) to WT (L) ###########

rawCountTable <- read.table("det_l_wt_l.bed")

cts <- as.matrix(rawCountTable)

########################"DESeq2 matrix #####################################################################


coldata <- read.csv("design_file_det_l_wt_l.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "WT")
sizeFactors(dds) <- c(1.00, 0.91,0.37, 0.32)
dds <- DESeq(dds, fitType='local')

res5 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)
table5 <- subset(res5, select=c("log2FoldChange", "padj", "pvalue"))
table5 <- as.data.frame(table5)
table5 <- na.omit(table5)
mcols(res5, use.names=TRUE)

ordered_res5 <- res5[order(res5$pvalue),]
summary(ordered_res5)


#############################################################################################
##########################Output tables with Wald test p-values##############################
#############################################################################################

write.table(as.data.frame(ordered_res1), col.names=NA, row.names = TRUE, sep= "\t", file="results_DE_WT_D_WT_L.tsv")
write.table(as.data.frame(ordered_res2), col.names=NA, row.names = TRUE, sep= "\t", file="results_DE_DET_D_DET_L.tsv")
write.table(as.data.frame(ordered_res3), sep= "\t", col.names=NA, row.names = TRUE, file="results_DE_DET_D_WT_D.tsv")
write.table(as.data.frame(ordered_res4), sep= "\t", col.names=NA, row.names = TRUE, file="results_DE_DET_D_WT_L.tsv")
write.table(as.data.frame(ordered_res5), sep= "\t", col.names=NA, row.names = TRUE, file="results_DE_DET_L_WT_L.tsv")




#########################################################################################################################
##########################MA plotsp-adj threshold 0.01 (1% false positives)##############################################
#########################################################################################################################

DESeq2::plotMA(res1, alpha = 0.01, ylim = c(-4, 4), main="MA plot WT Light : WT Dark", ylab= "log fold change WT Light : WT Dark")
DESeq2::plotMA(res2, alpha = 0.01, ylim = c(-4, 4), main="MA plot det1-1 Dark : det1-1 Light", ylab= "log fold change det1-1 Dark : det1-1 Light")
DESeq2::plotMA(res3, alpha = 0.01, ylim = c(-8, 6), main="MA plot det1-1 Dark : WT Dark", ylab= "log fold change det1-1 Dark : WT Dark")
DESeq2::plotMA(res4, alpha = 0.01, ylim = c(-8, 6), main="MA plot det1-1 Dark : WT Light", ylab= "log fold change det1-1 Dark : WT Light")
DESeq2::plotMA(res5, alpha = 0.01, ylim = c(-8, 6), main="MA plot det1-1 light : WT Light", ylab= "log fold change det1-1 light : WT Light")





