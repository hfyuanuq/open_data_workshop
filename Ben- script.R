Appendix 3.2 DE-Seq2 differential expression analysis R Script
##Combining the counts from both lanes of the CEL-Seq2 run
##For dataset 1(lane 1) make countsTable
countsTable <- read.delim(fle.choose(), header=TRUE, colClasses=c("character", rep("numeric", 36)))
#"numeric," followed by number refers to number of samples
rownames(countsTable) <- countsTable$Sample
countsTable <- countsTable[,-1] #says that the frst row is names not variables
##For dataset 2(lane 2) make countsTable2
countsTable2 <- read.delim(fle.choose(), header=TRUE, colClasses=c("character",rep("numeric", 36)))
rownames(countsTable2) <- countsTable2$Sample
countsTable2 <- countsTable2[,-1]
#To combine counts
temp <- cbind(countsTable, countsTable2)
countsall <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))
#To merge samples from di???erent matrices
challenge <- merge(as.data.frame(countsTable), as.data.frame(countsTable2), by="row.names", sort=FALSE)
##To print/write table
write.table(countsall, "c:/mydata.tab", sep="\t")
####################################################################################
#Load the packages into the R environment
library(DESeq2)
library(gplots)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
##Import fnal dataset(ERCC rows removed)
counts<-read.table(fle.choose(),header=TRUE,row.names=1,sep="\t")Deciphering the genomic toolkit underlying animal-bacteria interactions
182
##If gene names are in column 1 (not rownames) - convert to rownames
#counts_rename<-data.frame(counts[,-1],row.names=counts[,1])
#counts<-counts_rename
#rm(counts_rename)
#Use DESeq to prepare a data.frame of your full experimental design, this goes in as the colData input for
creating the DESeq2 object
treatment<-c(rep("foreign", 9),rep("healthy", 9),rep("control", 9),rep("unhealthy", 9)
             ,rep("foreign", 5),rep("healthy", 6),rep("control", 6),rep("unhealthy", 6))
time<-c(rep("0h", 3),rep("2h", 3),rep("8h", 3),rep("0h", 3),rep("2h", 3),rep("8h", 3),rep("0h", 3),rep("2h", 3),rep("8h",
                                                                                                                     3)
        ,rep("0h", 3),rep("2h", 3),rep("8h", 3),rep("0h", 2),rep("2h", 2),rep("8h", 1),rep("0h", 2),rep("2h", 2),rep("8h", 2)
        ,rep("0h", 2),rep("2h", 2),rep("8h", 2),rep("0h", 2),rep("2h", 2),rep("8h", 2))
design=data.frame(
  row.names=colnames(counts),
  condition=treatment,
  type=time)
##Create a DeSeq2 object, followed by ftting GLM to data
#flterstat: the flter statistic, here the average number of counts per gene across all samples, irrespective of
sample annoation,
#pvalue: the test p-values
dds <- DESeqDataSetFromMatrix(countData=counts, colData=design, design=~ type+condition) #creates a
DESeq2 object from counts matrix, experimental design and design forumula
colData(dds) #view the experimental design component of DESeq2 object, note that colData specifes information about the samples and experimental design. The frst column of colData would line up with the frst column
of counts matrix (i.e. genes)
as.data.frame(colData(dds)) #comprehensive view the experimental design component of DESeq2 object
design(dds) # view formula, condition, the variable of interest is at the end of the formula
##################Exploratory analysis and visualization#################
###Pre-fltering the dataset###
##Reduce the size of the object to increase the speed of the analyses by removing rows which little to no
information about gene expression.
nrow(dds) #lists number of rows183
Appendices
#overwrites existing DESeq2 object with an object in which genes with less than one read across all samples
(rowsum >1) are removed from the counts matrix
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
##Counts data are transformed to stabilise the variance across the mean. These become approximately
homoskedastic, and can be used directly for computing distances between samples and making PCA plots
#Sequencing depth correction is done automatically for the rlog and varianceStabilizingTransformation
#rlog and VST are meant for raw counts, use log2 transformation for normalised counts (after estimateSizeFactors)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
#export table
write.table(as.data.frame(assay(vsd)),fle='challenge-DESeq2-vst-transformed-counts.txt', sep='\t')
###Sample distances###
#Transpose the matrix of values using t, because the dist function expects the di???erent samples to be rows of
its argument, and genes to be columns.
sampleDists <- dist( t( assay(vsd) ) )
sampleDists
#Set up matrix
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$condition, vsd$type, sep="-" )
colnames(sampleDistMatrix) <- NULL
###PCA plot###
plotPCA(vsd, intgroup = c("condition", "type"))
library("ggplot2")
data <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  geom_text(aes(label=colnames(sampleDistMatrixVSD)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
#to change ggplot colour, use: scale_colour_brewer(palette="Set1") or change scale_colour_discrete(1 = 40)
# # # # # # # # # # # # # # # # # # # # # # # # # # # D i f f e r e n t i a l g e n e e x p r e s s i o nDeciphering the genomic toolkit underlying animal-bacteria interactions
184
analysis###############################################
dds$group <- factor(paste0(dds$type, dds$condition))#adds an additional column to colData that groups treatments for easier testing
design(dds) <- ~ group #Formula
dds <- DESeq(dds) #estimates size factors (control for di???erences in library size, dispersions, genewise dispersions, mean-dispersion relationship, fnal dispersion estimates, ftting GLM
resultsNames(dds)
## Plot of dispersion estimates (Diagnostic plot)
plotDispEsts(dds)
#export normalised data (DESeq2 normalisation).
write.table(counts(dds,normalized=TRUE),fle="challengeUMI_DESeq2normalized_counts.tab")
############The actual di???erential expression bit starts here##################
###First comparison with breakdown, subsequent comparisons will just be code######
## 2h, Foreign vs Control ##
# Fold change calculation: numerator followed by denominator
res2FC <- results(dds, contrast=c("group", "2hforeign", "2hcontrol"))
mcols(res2FC, use.names=TRUE)
# Merge fold change data with normalised counts
res2FCdata <- merge(as.data.frame(res2FC), as.data.frame(counts(dds, normalized=TRUE)), by="row.names",
                    sort=FALSE)
write.table(res2FCdata, fle="2h_ForeignVsControl_results.tab", sep = "\t")
# Genes with adjusted p-value less than 0.1 (10% false positive, padj=0.1)
sum(res2FC$padj < 0.1, na.rm=TRUE)
res2FCSig <- subset(res2FC, padj < 0.1)
# Sort it by the log2 fold change estimate to get the signifcant genes with the strongest upregulation
res2FCSig_ordered <- res2FCSig[ order( -res2FCSig$log2FoldChange, decreasing=TRUE), ]
write.table(res2FCSig_ordered, fle="2h_ForeignVsControl__padj0.1.tab", sep = "\t")
#### Subsequent comparisons begin here ####185
Appendices
### 2h ###
## 2h, Healthy vs Control ##
res2HC <- results(dds, contrast=c("group", "2hhealthy", "2hcontrol"))
sum(res2HC$padj < 0.1, na.rm=TRUE)
res2HCSig <- subset(res2HC, padj < 0.1)
res2HCSig_ordered <- res2HCSig[ order( -res2HCSig$log2FoldChange, decreasing=TRUE ), ]
write.table(res2HCSig_ordered, fle="2h_HealthyVsControl__padj0.1.tab", sep = "\t")
## 2h, Unhealthy vs Control ##
res2UC <- results(dds, contrast=c("group", "2hunhealthy", "2hcontrol"))
sum(res2UC$padj < 0.1, na.rm=TRUE)
res2UCSig <- subset(res2UC, padj < 0.1)
res2UCSig_ordered <- res2UCSig[ order( -res2UCSig$log2FoldChange, decreasing=TRUE ), ]
write.table(res2UCSig_ordered, fle="2h_unHealthyVsControl__padj0.1.tab", sep = "\t")
### 8h ###
## 8h, Healthy vs Control ##
res8HC <- results(dds, contrast=c("group", "8hhealthy", "8hcontrol"))
sum(res8HC$padj < 0.1, na.rm=TRUE)
res8HCSig <- subset(res8HC, padj < 0.1)
res8HCSig_ordered <- res8HCSig[ order( -res8HCSig$log2FoldChange, decreasing=TRUE ), ]
write.table(res8HCSig_ordered, fle="8h_HealthyVsControl__padj0.1.tab", sep = "\t")
## 8h, Unhealthy vs Control ##
res8UC <- results(dds, contrast=c("group", "8hunhealthy", "8hcontrol"))
sum(res8UC$padj < 0.1, na.rm=TRUE)
res8UCSig <- subset(res8UC, padj < 0.1)
res8UCSig_ordered <- res8UCSig[ order( -res8UCSig$log2FoldChange, decreasing=TRUE ), ]
write.table(res8UCSig_ordered, fle="8h_unHealthyVsControl__padj0.1.tab", sep = "\t")
##To view results without indepenedent fltering
#resNoFilt <- results(dds, contrast=c("group", "2hforeign", "2hcontrol"), independentFiltering=FALSE)
#addmargins(table(fltering=(res$padj < .1), noFiltering=(resNoFilt$padj < .1)))
#sum(resNoFilt$padj < 0.1, na.rm=TRUE)
### Diagnostic plots ###Deciphering the genomic toolkit underlying animal-bacteria interactions
186
## MA-plot The log2 fold change for a particular comparison is plotted on the y-axis and the average of the
counts normalized by size factor is shown on the x-axis.
##("M" for minus, because a log ratio is equal to log minus log, and "A" for average)
plotMA(res2FC, ylim=c(-5,5))
## Examine independent fltering
attr(res2FC, "flterThreshold")
plot(attr(res2FC,"flterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")