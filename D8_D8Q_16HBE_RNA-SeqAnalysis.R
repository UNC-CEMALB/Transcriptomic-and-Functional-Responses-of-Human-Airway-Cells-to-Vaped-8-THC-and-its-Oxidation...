
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
install.packages("pheatmap")
BiocManager::install("org.Hs.eg.db") #change this if using a species other than human
install.packages("ComplexHeatmap")
install.packages("tibble")
install.packages("dplyr")
BiocManager::install("fgsea")
install.packages("tidyverse")

library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(fgsea)
library(ggplot2)

#Loading in data and metadata
setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/CBD and D8 Project/RNA-Seq/CBD D8 RNA Seq/My Analysis")
cts <- read.csv('raw_counts_D8_D8Q_12.csv', header = TRUE, sep = ",")
coldata <- read.csv('D8_D8Q_12_coldata.csv', header = TRUE, sep = ",")

#run the following lines and check that order of samples is the same in cts and coldata
head(cts)
coldata

#constructing a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition, tidy = TRUE)
dds

#pre-filtering genes
smallestGroupSize <- 3 #change this to the minimum number of samples you had in any of your treatment groups
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize #this requires genes to have 3 or more samples with counts of 10 or more. 
dds <- dds[keep,]

#specifying the reference level (control)
dds$condition <- relevel(dds$condition, ref = "VC")

#differential expression analysis for all groups combined 
dds <- DESeq(dds)

#transformed values
vsd <- vst(dds, blind=FALSE) #variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #regularized log transformation
ntd <- normTransform(dds) #normalized counts transformation
head(assay(vsd), 3)

#PCA of all groups
tiff(filename = "D8THC_D8THCQ_PCA_12.tiff", width = 1000, height = 2000, res = 300)
plotPCA(vsd, intgroup=c("condition")) #run this line only to visualize in R without saving
dev.off()

library(ggrepel)
library(ggfortify)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- factor(pcaData$condition, levels=c("VC", "D8THC", "D8THCQ"), labels=c("VC", "D8THC", "D8THCQ"))

#better looking PCA
png(filename = "D8THC_D8THCq_12hr_PCA.png", width = 5500, height = 4000, res = 900)
pca <- ggplot(pcaData, aes(PC1, PC2))
pca + geom_point(size=3, aes(color=condition, fill=condition, group=condition)) + 
  ggtitle("Principle Component Analysis 12 hours") + 
  xlab(paste0("PC1: ", percentVar[1],"% varience")) + 
  ylab(paste0("PC2: ", percentVar[2],"% varience")) + 
  theme(legend.title=element_blank()) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("VC" = "lightblue", "D8THC" = "#B59FDC", "D8THCQ" = "#6D2E46")) +
  coord_fixed() 
dev.off()


#Getting human gene symbols
library("org.Hs.eg.db") #downloading human database, need to change if different species
columns(org.Hs.eg.db)

#Extracting results for each separate comparison
#comparison 1 (D8THC vs VC)
res1 <- results(dds, contrast=c("condition","D8THC","VC")) #reference or control should be second
resOrdered1 <- res1[order(res1$padj),] #ordering by adjusted p-value
res.df1 <- as.data.frame(resOrdered1)
res.df1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df1), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df1, file="D8THCvsVC_12.csv") #writing results file for this specific comparison
res05_1 <- results(dds, contrast=c("condition","D8THC","VC"), alpha=0.05)
summary(res05_1) #use this output to see an overview of results

#comparison 2 (D8THCQ vs VC)
res2 <- results(dds, contrast=c("condition","D8THCQ","VC")) #reference or control should be second
resOrdered2 <- res2[order(res2$padj),] #ordering by adjusted p-value
res.df2 <- as.data.frame(resOrdered2)
res.df2$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df2), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df2, file="D8THCQvsVC_12.csv") #writing results file for this specific comparison
res05_2 <- results(dds, contrast=c("condition","D8THCQ","VC"), alpha=0.05)
summary(res05_2) #use this output to see an overview of results

#comparison 3 (D8THCQ vs D8THC)
res3 <- results(dds, contrast=c("condition","D8THCQ","D8THC")) #reference or control should be second
resOrdered3 <- res3[order(res3$padj),] #ordering by adjusted p-value
res.df3 <- as.data.frame(resOrdered3)
res.df3$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df3), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df3, file="D8THCQvsD8THC_12.csv") #writing results file for this specific comparison
res05_3 <- results(dds, contrast=c("condition","D8THCQ","D8THC"), alpha=0.05)
summary(res05_3) #use this output to see an overview of results

#Making results file of only significant genes for each comparison, depending on results may need to alter these cut-offs
#Then making a vector of significant gene probe IDs to use in heatmap
res_sig1 <- res1[which(res1$padj < 0.05 & abs(res1$log2FoldChange) >= 1), ] 
sig_D8THCvsVC <- rownames(res_sig1)
res_sig1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_sig1), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig1, file="D8THCvsVC_sig_12.csv")
summary(res_sig1)

res_sig2 <- res2[which(res2$padj < 0.05 & abs(res2$log2FoldChange) >= 1), ] 
sig_D8THCQvsVC <- rownames(res_sig2)
res_sig2$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_sig2), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig2, file="D8THCQvsVC_sig_12.csv")
summary(res_sig2)

res_sig3 <- res3[which(res3$padj < 0.05 & abs(res3$log2FoldChange) >= 1), ] 
sig_D8THCQvsD8THC <- rownames(res_sig3)
res_sig3$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_sig3), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig3, file="D8THCQvsD8THC_sig_12.csv")
summary(res_sig3)

#Combining vectors of significant genes for each comparison and then removing duplicate probe IDs
all_sig_genes <- c(sig_D8THCvsVC, sig_D8THCQvsVC, sig_D8THCQvsD8THC)
all_sig_genes_unique <- unique(all_sig_genes)

#Making matrix of normalized counts to use in heatmap, normalized counts are typically used for this type of visualization
normcounts <- counts(dds, normalized=TRUE)
head(normcounts)

#filtering normalized counts by vector of significant genes previously made
norm_sigonly <- subset(normcounts, rownames(normcounts) %in% all_sig_genes_unique)

#Setting up to make heatmap
hmtable <- norm_sigonly
head(hmtable)
dim(hmtable)

#scaling data to obtain z-scores
hmtable_transposed <- t(hmtable)
head(hmtable_transposed)
hmtable_scaled <- scale(hmtable_transposed)
hmtable <- t(hmtable_scaled)

#plotting and saving heatmap 
png(filename = "D8_12_HM.png", width = 7500, height = 6000, res = 900)
Heatmap(hmtable,
        #column_dend_reorder = FALSE,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(Group = coldata$condition, simple_anno_size = unit(.4, "cm"),
                                           show_annotation_name = FALSE,                                  
                                           col = list(Group = c("VC" = "lightblue", "D8THC" = "#B59FDC", "D8THCQ" = "#6D2E46")),
                                           annotation_legend_param = list(Group = list(labels = c("VC", "D8THC", "D8THCQ"), at = c("VC", "D8THC", "D8THCQ")))))
dev.off()

#Making Venn diagram
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library("ggVennDiagram")

#ggVennDiagram works for up to 7 comparisons, just add more comparisons to the list function below
#using the vectors of significant genes made previously
venn_list <- list(#A = sig_D8THCvsVC, 
                  B = sig_D8THCQvsVC, 
                  C = sig_D8THCQvsD8THC)

tiff(filename = "D8THC_D8THCQ_Venn.tiff", width = 2000, height = 1000, res = 300)
ggVennDiagram(venn_list, category.names = c("D8THCQ vs. VC", "D8THCQ vs. D8THC")) +
  scale_x_continuous(expand = expansion(mult = .1)) #expanding x axis so labels aren't cut off
dev.off()

library(gplots)
venn(venn_list)
my_venn <- venn(venn_list)
class(my_venn)
library(pryr)
names(attributes(my_venn))
venn_intersections <- attr(x = my_venn, "intersections")
vennDF <- print(as.data.frame(do.call(cbind, venn_intersections)))
write.csv(vennDF, file ="EulerIntersections_VapedCBD.csv")

library(eulerr)

fit1 <- euler(venn_list)
fit1
fit2 <- euler(venn_list, shape = "ellipse")
fit2
plot(fit2)

custom_colors <- c("D8THCQ vs. VC" = "#6D2E46",   
                   "D8THCQ vs. D8THC" = "grey")

png(filename = "D8_12_Sig_Euler.png", width = 3500, height = 3500, res = 900)
plot(fit2, 
     legend = list(side = "bottom"),
     fills = custom_colors,
     quantities = list(type = c("counts")))
dev.off()

#to get gene lists from certain intersections in your Venn diagram
intersect <- intersect(sig_D8THCQvsVC, sig_D8THCQvsD8THC) #getting overlapping probe IDs
intersect <- as.data.frame(intersect) #making it a dataframe
colnames(intersect) = c("GeneID")
intersect$GeneName <- mapIds(org.Hs.eg.db, keys = intersect$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
intersect$GeneName

#see the following for more customization of the ggVennDiagram
  #https://cran.r-project.org/web/packages/ggVennDiagram/readme/README.html
  #https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html#:~:text=If%20you%20use%20long%20category,trick%20to%20expand%20x%20axis.

#_________VOLCANO PLOTS_________
library(ggrepel)

# Ensure your dataframe has the correct structure
res.df2$gene <- ifelse(is.na(res.df2$symbol) | res.df2$symbol == "", rownames(res.df2), res.df2$symbol)

# Define threshold for significance
log2FC_threshold <- 1  # Adjust as needed
pval_threshold <- 0.05

# Add a column for significance
res.df2$Significance <- "Not Significant"
res.df2$Significance[res.df2$log2FoldChange > log2FC_threshold & res.df2$padj < pval_threshold] <- "Upregulated"
res.df2$Significance[res.df2$log2FoldChange < -log2FC_threshold & res.df2$padj < pval_threshold] <- "Downregulated"

# Convert to factor for proper color mapping
res.df2$Significance <- factor(res.df2$Significance, levels = c("Upregulated", "Downregulated", "Not Significant"))

# Create the volcano plot
volcano_plot <- ggplot(res.df2, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  xlab("Log2 Fold Change") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot of DEGs") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

# Add labels for top genes
top_genes <- res.df2[order(res.df2$padj), ][1:10, ]  # Adjust number as needed
volcano_plot <- volcano_plot + 
  geom_label_repel(
    data = top_genes, aes(label = gene), 
    size = 3, 
    box.padding = 0.5, 
    max.overlaps = 10, 
    fill = "white",  # White background for better contrast
    color = "black", # Black text for readability
    alpha = 0.75,    # Slight transparency to blend naturally
    segment.color = "gray50" # Connecting line to the point
  )

# Print plot

print(volcano_plot)
ggsave("D8QvsVC_volcano_24.png", plot = volcano_plot, width = 853/96, height = 604/96, dpi = 900, units = "in")
