library("DESeq2")
library("ggplot2")
library("gplots")
library("dplyr")
library("ggrepel")
library("genefilter")
library("RColorBrewer")
library("GenomicAlignments")
library("IHW")
library("BiocParallel")
library("data.table")
library("limma")
library("geneplotter")
register(MulticoreParam(4))
library("amap")

################################################ DE analysis based on HTSeq counts for det1 data project ###########################

#### Read the HTseq counts files 

basedir <- "/export/home1/users/mpb/amohamed/Mounted_dir/Project_mpb_pourlemomentnon/barneche/Ouardia/chip_seq_data/det1_project/rna_seq_data/htseq_count_files"
setwd(basedir)
cntdir <- paste(basedir, "dark", sep="/")
pat <- ".tsv"
tophat.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# we choose the 'all' series
myfiles <- tophat.all
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*)_all_counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# inspect
head(data)
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# inspect and look at the top row names!
head(data)

# create a table like organized in the colData file
rawCountTable <- data

colnames(rawCountTable) <- c("det1-1_D1", "det1-1_D2", 
                             "det1-1_L1", "det1-1_L2", 
                             "WT_D1", "WT_D2",
                             "WT_L1", "WT_L2")

cts <- as.matrix(rawCountTable)

#### Read the design file #############

coldata <- read.csv("design_file_det1_rna_seq_dark.csv")
coldata <- as.data.frame(coldata)
coldata$replicate <- factor(coldata$replicate)


############## DDS object ##################

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)

dds$group <- factor(paste0(dds$treatment, dds$genotype))
design(dds) <- ~ group
dds <- DESeq(dds, minReplicatesForReplace = Inf)


####Box plot of Data ####
par(mar=c(9,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, col="lightblue")

####Model distribution plot of Data ####
plotDispEsts(dds)

############## Regularized log transformation of Data ##################
rld <- rlog(dds)

############## Checking the quality of Data ##################

####Spearman correlation Hierarchical clustering ####
plot(hcluster(dist(t(counts(dds))), method= "spearman", link="average"), cex=0.5, cex.main=1, main= "Hierarchical Clustering of Treated det1-1 Count Data")


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


colours = colorRampPalette(rev(brewer.pal(9, "Reds")))(256)
heatmap.2(sampleDistMatrix, trace="none", col=colours, margins = c(9, 9)) 

####Gene expression Heatmap plot####

topVarGenest <- head(order(rowVars(assay(rld)), decreasing=TRUE), 22319)
heatmap.2(assay(rld)[topVarGenest, ], scale="row",labRow = FALSE,
          trace="none", dendrogram="column",
          col = colorRampPalette(c("green","black","red"))(256), margins = c(10, 10), 
          main="Sample Distance Matrix")


######################################################################
##########################Contrasts###################################
######################################################################


############## Contrast det1-1 (D) vs det1-1 (L) ###########

res1 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE, contrast=c("group", "Ddet1-1", "Ldet1-1"))
table1 <- subset(res1, select=c("log2FoldChange", "padj", "pvalue"))
table1 <- as.data.frame(table1)
table1 <- na.omit(table1)
mcols(res1, use.names=TRUE)


ordered_res1 <- res1[order(res1$pvalue),]


############## Contrast WT (D) vs WT (L) ###########

res2 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE, contrast=c("group", "DWT", "LWT"))
table2 <- subset(res2, select=c("log2FoldChange", "padj", "pvalue"))
table2 <- as.data.frame(table2)
table2 <- na.omit(table2)
mcols(res2, use.names=TRUE)

ordered_res2 <- res2[order(res2$pvalue),]


############## Contrast det1-1 (D) vs WT (D) ###########

res3 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE, contrast=c("group", "Ddet1-1", "DWT"))
table3 <- subset(res3, select=c("log2FoldChange", "padj", "pvalue"))
table3 <- as.data.frame(table3)
table3 <- na.omit(table3)
mcols(res3, use.names=TRUE)

ordered_res3 <- res3[order(res3$pvalue),]


############## Contrast det1-1 (D) vs WT (L) ###########

res4 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE, contrast=c("group", "Ddet1-1", "LWT"))
table4 <- subset(res4, select=c("log2FoldChange", "padj", "pvalue"))
table4 <- as.data.frame(table4)
table4 <- na.omit(table4)
mcols(res4, use.names=TRUE)

ordered_res4 <- res4[order(res4$pvalue),]


############## Contrast det1-1 (L) vs WT (L) ###########

res5 <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE, contrast=c("group", "Ldet1-1", "LWT"))
table5 <- subset(res5, select=c("log2FoldChange", "padj", "pvalue"))
table5 <- as.data.frame(table5)
table5 <- na.omit(table5)
mcols(res5, use.names=TRUE)

ordered_res5 <- res5[order(res5$pvalue),]


#############################################################################################
##########################Output tables with Wald test p-values##############################
#############################################################################################

write.table(as.data.frame(ordered_res1), col.names=NA, row.names = TRUE, sep= "\t", file="res_DE_rna_seq_det_d_det_l.tsv")
write.table(as.data.frame(ordered_res2), col.names=NA, row.names = TRUE, sep= "\t", file="res_DE_rna_seq_wt_d_wt_l.tsv")
write.table(as.data.frame(ordered_res3), col.names=NA, row.names = TRUE, sep= "\t", file="res_DE_rna_seq_det_d_wt_d.tsv")
write.table(as.data.frame(ordered_res4), col.names=NA, row.names = TRUE, sep= "\t", file="res_DE_rna_seq_det_d_wt_l.tsv")
write.table(as.data.frame(ordered_res5), col.names=NA, row.names = TRUE, sep= "\t", file="res_DE_rna_seq_det_l_wt_l.tsv")


#########################################################################################################################
##########################MA plotsp-adj threshold 0.01 (1% false positives)####################################
#########################################################################################################################

DESeq2::plotMA(res1, alpha = 0.01, ylim = c(-14, 14), main="MA plot det1-1 Dark:det1-1 Light", ylab= "log fold change det1-1 Dark:det1-1 Light")
DESeq2::plotMA(res2, alpha = 0.01, ylim = c(-14, 14), main="MA plot WT Dark:WT Light", ylab= "log fold change WT Dark:WT Light")
DESeq2::plotMA(res3, alpha = 0.01, ylim = c(-14, 14), main="MA plot det1-1 Dark:WT Dark", ylab= "log fold change det1-1 Dark:WT Dark")
DESeq2::plotMA(res4, alpha = 0.01, ylim = c(-14, 14), main="MA plot det1-1 Dark:WT Light", ylab= "log fold change det1-1 Dark:WT Light")
DESeq2::plotMA(res5, alpha = 0.01, ylim = c(-14, 14), main="MA plot det1-1 Light:WT Light", ylab= "log fold change det1-1 Light:WT Light")




