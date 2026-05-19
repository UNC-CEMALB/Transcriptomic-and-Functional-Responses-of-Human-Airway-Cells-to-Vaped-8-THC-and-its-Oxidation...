# ==============================================================================
# SECTION: Package and Environment Setup
# ==============================================================================

# 1. Define required packages
required_packages <- c(
  "tidyverse", "dplyr", "stringr", "ggplot2", "purrr", "readxl", "data.table", 
  "psych", "corrplot", "openxlsx", "EnhancedVolcano", "eulerr", "ggVennDiagram", 
  "pheatmap", "RColorBrewer", "edgeR", "DESeq2", "AnnotationDbi", "org.Hs.eg.db", 
  "clusterProfiler", "enrichplot", "lobstr", "gplots"
)

# 2. Quietly handle installations
if (!suppressPackageStartupMessages(requireNamespace("BiocManager", quietly = TRUE))) {
  install.packages("BiocManager", quietly = TRUE)
}

# 3. Suppress output and resolve dependencies/missing packages
# This will find and install lobstr (from CRAN) and gplots automatically
suppressPackageStartupMessages({
  BiocManager::install(required_packages, update = TRUE, ask = FALSE, 
                       dependencies = TRUE, quietly = TRUE)
})

# 4. Silent, efficient loading
silence_output <- lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
})

# 5. Organism-specific setup
organism <- "org.Hs.eg.db"
suppressPackageStartupMessages(library(organism, character.only = TRUE, quietly = TRUE))

# 6. Final Clean-up: Remove setup variables to keep environment pristine
rm(required_packages, silence_output)

# --- Load required packages ---

# Tidyverse components
library(tidyverse)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(purrr)
library(readxl)
library(data.table)
library(openxlsx)

# Specialized statistical and visualization packages
library(psych)
library(corrplot)
library(EnhancedVolcano)
library(eulerr)
library(ggVennDiagram)
library(pheatmap)
library(RColorBrewer)

# Bioconductor and analysis packages
library(edgeR)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# Organism-specific database setup (Optional but good practice)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# Confirmation message
print("All required packages are now loaded.")


setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Jaspers Lab/CBD and D8 Project/RNA-Seq/Vaped D8 RNA-Seq/My Analysis")

cts <- read.csv('AllSamples_RawCounts_VapedD8RNASeq.csv', header = TRUE, sep = ",")#raw data # reading in RNAseq data
coldata <- read.csv('Metadata_AllSamples_RawCounts_VapedD8RNASeq.csv', header = TRUE, sep = ",")

counts <- cts

d0 <- DGEList(counts)# Create DGEList object
d0 <- calcNormFactors(d0)
cutoff <- 4 
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

snames <- colnames(counts) # Sample names
snames <- snames[2:length(snames)]

#removing replicate suffix for tech reps and making replicate vector
base_names <- gsub("_[0-9]+$", "", snames)
replicate <- factor(base_names)

#extracting product and donor info
product <- sapply(snames, function(x) {
  strsplit(x, "_")[[1]][3]  # Split by "_" and get the 3rd part
}) 
print(product)

group <- as.factor(product)
group

donor <- sapply(snames, function(x) {
  strsplit(x, "_")[[1]][1]  # Split by "_" and get the 1st part
})
print(donor)

batch <- as.factor(donor)
batch

plotMDS(d, col = as.numeric(group))

mm <- model.matrix(~0+group+batch)

y <- voom(d, mm, plot = T)

#performing duplicateCorrelation to account for tech reps
corfit <- duplicateCorrelation(y, mm, block = replicate)
fit <- lmFit(y, mm, block = replicate, correlation = corfit$consensus.correlation)
fit <- eBayes(fit)
top.table <- topTable(fit)
df_fit <- as.data.frame(fit$coefficients)


df <- as.data.frame(y$E)
adjusted_values <- fitted(fit)
df_adjusted <- as.data.frame(adjusted_values)
colnames(df_adjusted) <- colnames(y$E)
head(df_adjusted)
df_adjusted_no_2 <- df_adjusted[, !grepl("_2$", colnames(df_adjusted))]
head(df_adjusted_no_2)
df <- df_adjusted_no_2


## PCA/MDS with donor correction for all groups

row.names(df) <- y$genes$Geneid

Labels <- substr(colnames(df), 1, 18)
Product <- substr(colnames(df), 18, 18)
donor <- substr(colnames(df), 1, 8) # for batch correction

dfB <- limma::removeBatchEffect(df, donor)

group_colors <- c(rep("#8DD3C7", 4), rep("#FFFFB3", 4), rep("#BEBADA", 4), 
                  rep("#FB8072", 4), rep("#B3DE69", 4), rep("#FDB462", 4), rep("lightblue", 4))
group_labels <- c("Palmetto Distillate", "NYSW Distillate", "Palmetto Disposable", "Stiizy Disposable", "ST Strawberry Juice", "ST Green Apple Juice", "PG/VG Control")

# Define the shape codes (you can adjust the shape codes as needed)
shape_codes <- c(15, 16, 17, 18)  # Define 4 shapes, e.g., square, circle, triangle, diamond

# Repeat the shape codes in blocks of 4 shapes for every 4th sample
n_samples <- length(colnames(df)) 
group_shapes <- rep(shape_codes, length.out = n_samples)

# Ensure every 4th sample has the same shape (adjust by position)
group_shapes[seq(1, n_samples, by = 4)] <- shape_codes[1]  # Shape for 1, 5, 9, ...
group_shapes[seq(2, n_samples, by = 4)] <- shape_codes[2]  # Shape for 2, 6, 10, ...
group_shapes[seq(3, n_samples, by = 4)] <- shape_codes[3]  # Shape for 3, 7, 11, ...
group_shapes[seq(4, n_samples, by = 4)] <- shape_codes[4]  # Shape for 4, 8, 12, ...

#ensure order of samples correctly matches order of the colors and shapes
print(colnames(df)) #sample order
print(group_colors)
print(group_shapes)

plotMDS(df, pch = 10, cex = 1.2, labels = Labels, gene.selection = "pairwise", main="Original", col=group_colors)
mds_original <- plotMDS(df, pch = 10, cex = 1.2, labels = Labels, plot = FALSE, gene.selection = "pairwise", main="Original")

plotMDS(dfB, pch = 10, cex = 1.2, labels = Labels, gene.selection = "pairwise", main="After Donor Correction", col=group_colors)
mds_corrected <- plotMDS(dfB, pch = 10, cex = 1.2, labels = Labels, plot = FALSE, gene.selection = "pairwise", main="After Donor Correction")

x_original <- jitter(mds_original$x, amount = 0.05)  # Adjust 'amount' to control jittering
y_original <- jitter(mds_original$y, amount = 0.03)
x_corrected <- jitter(mds_corrected$x, amount = 0.05)
y_corrected <- jitter(mds_corrected$y, amount = 0.03)

# Convert MDS coordinates into a data frame for clustering
mds_original_df <- data.frame(Leading_logFC1 = x_original, Leading_logFC2 = y_original)
mds_corrected_df <- data.frame(Leading_logFC1 = x_corrected, Leading_logFC2 = y_corrected)

# Ensure product order is preserved
mds_original_df$Product <- factor(rep(group_labels, each = 4), levels = group_labels)  
mds_corrected_df$Product <- factor(rep(group_labels, each = 4), levels = group_labels)

# Ensure donor order is correct
mds_original_df$Donor <- factor(rep(1:4, length.out = nrow(mds_original_df)))
mds_corrected_df$Donor <- factor(rep(1:4, length.out = nrow(mds_corrected_df)))

# Run k-means clustering (adjust centers based on observed clusters)
set.seed(123)  # Ensures reproducibility
num_clusters <- 3  # Adjust based on observed clustering patterns
mds_original_df$Cluster <- factor(kmeans(mds_original_df[, c("Leading_logFC1", "Leading_logFC2")], 
                                         centers = num_clusters, nstart = 25)$cluster)
mds_corrected_df$Cluster <- factor(kmeans(mds_corrected_df[, c("Leading_logFC1", "Leading_logFC2")], 
                                          centers = num_clusters, nstart = 25)$cluster)

# Plot MDS with Data-Driven Ellipses (Before Batch Correction)
png(filename = "MDS_original_jitter_VapedD8.png", res = 900)
ggplot(mds_original_df, aes(x = Leading_logFC1, y = Leading_logFC2, color = Product, shape = Donor)) +
  geom_point(size = 3) +  
  stat_ellipse(aes(group = Cluster), level = 0.95, alpha = 0.2, color = "black") +  
  theme_bw() +  # Simple, clean theme
  labs(title = "MDS (Leading LogFC) with Data-Driven Clusters",
       x = "Leading logFC Dimension 1 (42%)",  
       y = "Leading logFC Dimension 2 (24%)") +  
  scale_color_manual(values = unique(group_colors)) +
  scale_shape_manual(values = unique(group_shapes)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) 
ggsave("MDS_original_VapedD8.png", 
       width = 650 / 100,  # Convert pixels to inches (900 DPI)
       height = 450 / 100, 
       dpi = 900)

# Plot MDS with Data-Driven Ellipses (After Donor Correction)
ggplot(mds_corrected_df, aes(x = Leading_logFC1, y = Leading_logFC2, color = Product, shape = Donor)) +
  geom_point(size = 3) +  
  stat_ellipse(aes(group = Cluster), level = 0.95, alpha = 0.2, color = "black") +  
  theme_bw() +  # Simple, clean theme
  labs(title = "MDS After Donor Correction with Data-Driven Clusters",
       x = "Leading logFC Dimension 1 (86%)",
       y = "Leading logFC Dimension 2 (5%)") +
  scale_color_manual(values = unique(group_colors)) +
  scale_shape_manual(values = unique(group_shapes)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) 
ggsave("MDS_corrected_VapedD8.png", 
       width = 650 / 100,  # Convert pixels to inches (900 DPI)
       height = 450 / 100, 
       dpi = 900)

##___________older plots________#
png(filename = "PCA_original_jitter.png", width = 8200, height = 5546, res = 900)
par(xpd = TRUE, mar = c(5, 6, 4, 10))
plot(x_original, y_original, 
     pch = group_shapes, 
     col = group_colors, 
     cex = 1.2,
     main = "Original with Jitter",
     cex.lab = 1.5, 
     cex.axis = 1.5)
legend("topright", 
       legend = group_labels,  # Product groups
       col = unique(group_colors), 
       pch = 19,  # Solid circles for color legend
       cex = 0.8,  
       title = "Product Type",
       xpd = TRUE,  # Allows plotting outside the figure region
       inset = c(-0.285, 0))
legend("bottomright", 
       legend = c("Donor 1", "Donor 2", "Donor 3", "Donor 4"),  # Donor groups
       col = "black",  # Shapes in black for clarity
       pch = unique(group_shapes),  # Uses the shape codes assigned to donors
       cex = 0.8,
       title = "Donor",
       xpd = TRUE,  # Allows plotting outside the figure region
       inset = c(-0.145, 0.425))
dev.off()

png(filename = "MDS_donorcorrected_jitter.png", width = 8200, height = 5546, res = 900)
par(xpd = TRUE, mar = c(5, 6, 4, 10))
plot(x_corrected, y_corrected, 
     pch = group_shapes, 
     col = group_colors, 
     cex = 1.2, 
     #xlab = x_label_after,
     #ylab = y_label_after, 
     main = "Donor-Corrected with Jitter",
     cex.lab = 1.5, 
     cex.axis = 1.5) 
legend("topright", 
       legend = group_labels,  # Product groups
       col = unique(group_colors), 
       pch = 19,  # Solid circles for color legend
       cex = 0.8,  
       title = "Product Type",
       xpd = TRUE,  # Allows plotting outside the figure region
       inset = c(-0.285, 0))
legend("bottomright", 
       legend = c("Donor 1", "Donor 2", "Donor 3", "Donor 4"),  # Donor groups
       col = "black",  # Shapes in black for clarity
       pch = unique(group_shapes),  # Uses the shape codes assigned to donors
       cex = 0.8,
       title = "Donor",
       xpd = TRUE,  # Allows plotting outside the figure region
       inset = c(-0.145, 0.425))
dev.off()

#_______ plotting actual PCA not MDS_______#

pca_df <- data.frame(PC1 = pca_result$x[,1], 
                     PC2 = pca_result$x[,2], 
                     Product = factor(rep(group_labels, each = 4)))  # Ensure grouping

# Determine clusters based on PCA structure
set.seed(123)  # Ensures reproducibility
num_clusters <- 3  # Adjust based on observed clustering
kmeans_result <- kmeans(pca_result$x[, 1:2], centers = num_clusters, nstart = 25)

# Convert PCA results into a data frame
pca_df <- data.frame(PC1 = pca_result$x[,1], 
                     PC2 = pca_result$x[,2], 
                     Product = factor(rep(group_labels, each = 4)),  # Product groups for color
                     Donor = factor(rep(1:4, length.out = nrow(pca_result$x))),  # Donors for shape
                     Cluster = as.factor(kmeans_result$cluster))  # Use detected clusters

# Plot PCA with:
# - Product colors
# - Donor shapes
# - Ellipses around clusters
ggplot(pca_df, aes(x = PC1, y = PC2, color = Product, shape = Donor)) +
  geom_point(size = 3) +  # Points
  stat_ellipse(aes(group = Cluster, fill = Cluster), level = 0.95, alpha = 0.2, color = "black") +  # Cluster ellipses
  theme_minimal() +
  labs(title = "PCA with Data-Driven Clusters",
       x = paste0("PC1 (", pca_variance_percent[1], "%)"),
       y = paste0("PC2 (", pca_variance_percent[2], "%)")) +
  scale_color_manual(values = unique(group_colors)) +  # Keep original product colors
  scale_shape_manual(values = unique(group_shapes))  # Keep original donor shapes

# Apply jitter directly to PCA values to prevent overlapping points
set.seed(123)  # Ensure reproducibility
num_clusters <- 3  # Adjust based on PCA clustering pattern
kmeans_corrected <- kmeans(pca_result_after$x[, 1:2], centers = num_clusters, nstart = 25)

# Convert donor-corrected PCA results into a data frame with jitter
pca_corrected_df <- data.frame(
  PC1 = jitter(pca_result_after$x[,1], amount = 8),  # Jitter X-axis
  PC2 = jitter(pca_result_after$x[,2], amount = 6),  # Jitter Y-axis
  Product = factor(rep(group_labels, each = 4)),  # Product types for color
  Donor = factor(rep(1:4, length.out = nrow(pca_result_after$x))),  # Donors for shape
  Cluster = as.factor(kmeans_corrected$cluster)  # Clusters based on PCA
)

pca_corrected_df$Product <- factor(pca_corrected_df$Product, levels = unique(group_labels)) 
color_map <- setNames(group_colors, unique(group_labels))  # Ensures correct color mapping

# Plot PCA with jittered values, product colors, donor shapes, and ellipses
ggplot(pca_corrected_df, aes(x = PC1, y = PC2, color = Product, shape = Donor)) +
  geom_point(size = 3) +  # PCA points with donor shapes
  stat_ellipse(aes(group = Cluster), level = 0.95, alpha = 0.2, color = "black") + 
  geom_text_repel(aes(label = rownames(pca_corrected_df)), size = 3, max.overlaps = 30) +  # Sample labels# Ellipses around clusters
  theme_minimal() +
  labs(
    title = "Donor-Corrected PCA with Data-Driven Clusters (Jitter Applied)",
    x = paste0("PC1 (", pca_variance_percent_after[1], "%)"),
    y = paste0("PC2 (", pca_variance_percent_after[2], "%)")
  ) +
  scale_color_manual(values = unique(group_colors)) +  # Keep original product colors
  scale_shape_manual(values = unique(group_shapes))  # Keep original donor shapes

#---------------------------------------------------------#
#Making comparisons or Distillate AND disposable, Juice, and PGVG Control
type <- sapply(snames, function(x) {
  strsplit(x, "_")[[1]][4]  # Split by "_" and get the 3rd part
}) 
print(type)

type <- as.factor(type)
type

base_names <- gsub("_[0-9]+$", "", snames)
replicate <- factor(base_names)

mm <- model.matrix(~0+type+batch)

y <- voom(d, mm, plot = T)

#performing duplicateCorrelation to account for tech reps
corfit <- duplicateCorrelation(y, mm, block = replicate)
fit <- lmFit(y, mm, block = replicate, correlation = corfit$consensus.correlation)

contrast.matrix <- makeContrasts(
  Dis_vs_Con = typeDis - typeCon,
  Ju_vs_Con = typeJu - typeCon,
  Dis_vs_Ju = typeDis - typeJu,
  levels = colnames(mm)
)

fit_contrast <- contrasts.fit(fit, contrast.matrix)
fit_contrast <- eBayes(fit_contrast)

#Extracting each comparison
top.table_Dis_vs_Con <- topTable(fit_contrast, coef = "Dis_vs_Con", adjust.method = "BH", number = Inf)
top.table_Ju_vs_Con <- topTable(fit_contrast, coef = "Ju_vs_Con", adjust.method = "BH", number = Inf)
top.table_Dis_vs_Ju <- topTable(fit_contrast, coef = "Dis_vs_Ju", adjust.method = "BH", number = Inf)

# Save or view results
write.csv(top.table_Dis_vs_Con, "Dis_vs_Con.csv")
write.csv(top.table_Ju_vs_Con, "Ju_vs_Con.csv")
write.csv(top.table_Dis_vs_Ju, "Dis_vs_Ju.csv")

#Making results file of only significant genes for each comparison
#Then making a vector of significant gene probe IDs to use in heatmap
res1 <- top.table_Dis_vs_Con
res1$symbol <- mapIds(org.Hs.eg.db, keys = res1$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
res2 <- top.table_Ju_vs_Con
res2$symbol <- mapIds(org.Hs.eg.db, keys = res2$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
res3 <- top.table_Dis_vs_Ju
res3$symbol <- mapIds(org.Hs.eg.db, keys = res3$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names

res_sig1 <- res1[which(res1$adj.P.Val < 0.05 & abs(res1$logFC) >= 1), ] 
summary(res_sig1)
sum(res_sig1$logFC > 0, na.rm = TRUE)
sum(res_sig1$logFC < 0, na.rm = TRUE)
write.csv(res_sig1, file="Dis_vs_Con_sig.csv")
sig_Dis_vs_Con <- res_sig1$ID  

res_sig2 <- res2[which(res2$adj.P.Val < 0.05 & abs(res2$logFC) >= 1), ] 
summary(res_sig2)
sum(res_sig2$logFC > 0, na.rm = TRUE)
sum(res_sig2$logFC < 0, na.rm = TRUE)
write.csv(res_sig2, file="Ju_vs_Con_sig.csv")
sig_Ju_vs_Con <- res_sig2$ID 

res_sig3 <- res3[which(res3$adj.P.Val < 0.05 & abs(res3$logFC) >= 1), ] 
summary(res_sig3)
sum(res_sig3$logFC > 0, na.rm = TRUE)
sum(res_sig3$logFC < 0, na.rm = TRUE)
write.csv(res_sig3, file="Dis_vs_Ju_sig.csv")
sig_Dis_vs_Ju <- res_sig3$ID 


all_sig_genes <- c(sig_Dis_vs_Con, sig_Ju_vs_Con, sig_Dis_vs_Ju)
all_sig_genes_unique <- unique(all_sig_genes)

#Making matrix of normalized counts to use in heatmap
normalized_counts <- as.data.frame(y$E)       # Extract normalized counts
rownames(normalized_counts) <- y$genes$ID

base_names <- gsub("_[0-9]+$", "", colnames(normalized_counts))
unique_samples <- unique(base_names)

#average normalized counts for technical replicates
averaged_counts <- matrix(0, nrow = nrow(normalized_counts), ncol = length(unique_samples))
rownames(averaged_counts) <- rownames(normalized_counts)
colnames(averaged_counts) <- unique_samples

for (sample in unique_samples) {
  replicate_cols <- which(base_names == sample)  # Columns corresponding to this sample
  averaged_counts[, sample] <- rowMeans(normalized_counts[, replicate_cols, drop = FALSE])
}
head(averaged_counts)

normcounts <- averaged_counts

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
donors <- c("DD001Rp1", "DD004Op1", "DD008Rp1", "DD021Rp1")
coldata <- read.csv('Metadata_notechrep_RawCounts_VapedD8RNASeq.csv', header = TRUE, sep = ",")

library(ComplexHeatmap)

png(filename = "Sig_Genes_AllGroups_HM.png", width = 6000, height = 6000, res = 900)
Heatmap(hmtable,
        #column_dend_reorder = FALSE,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(
          Product = coldata$Product, 
          Donor = coldata$Donor,
          Type = coldata$Product_type,
          simple_anno_size = unit(.4, "cm"),
          show_annotation_name = FALSE,  
          col = list(Product = c("A" = "#8DD3C7", "C" = "#FFFFB3", "D" = "#BEBADA", "E" = "#FB8072", "F" = "#B3DE69", "G" = "#FDB462", "H" = "#B59FDC"), 
                     Donor = c("DD001Rp1" = "#F0027F","DD004Op1" = "dodgerblue","DD008Rp1" = "#FDDA0D","DD021Rp1" = "#00CED1"),
                     Type = c("Distillate" = "#66C2A5", "Disposable" = "#FC8D62", "Juice" = "#8DA0CB", "PGVG" = "lightblue")),
          annotation_legend_param = list(Type = list(labels = c("Distillate", "Disposable", "Juice","PGVG"), at = c("Distillate", "Disposable", "Juice","PGVG")))))
dev.off()


library(magrittr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Create and store heatmap with fixed clustering
HM <- Heatmap(hmtable,
              km = 6,  # Ensures 6 clusters
              heatmap_legend_param = list(title = "Z-score"),
              column_names_gp = gpar(fontsize = 11),
              row_names_gp = gpar(fontsize = 6),
              show_row_names = FALSE,
              top_annotation = HeatmapAnnotation(
                Product = coldata$Product, 
                Donor = coldata$Donor,
                Type = coldata$Product_type,
                simple_anno_size = unit(.4, "cm"),
                show_annotation_name = FALSE,  
                col = list(Product = c("A" = "#8DD3C7", "C" = "#FFFFB3", "D" = "#BEBADA", "E" = "#FB8072", "F" = "#B3DE69", "G" = "#FDB462", "H" = "#B59FDC"), 
                           Donor = c("DD001Rp1" = "#F0027F","DD004Op1" = "dodgerblue","DD008Rp1" = "#FDDA0D","DD021Rp1" = "#00CED1"),
                           Type = c("Distillate" = "#66C2A5", "Disposable" = "#FC8D62", "Juice" = "#8DA0CB", "PGVG" = "lightblue")),
                annotation_legend_param = list(Type = list(labels = c("Distillate", "Disposable", "Juice","PGVG"), at = c("Distillate", "Disposable", "Juice","PGVG")))))

# **Force heatmap to be drawn and finalized**
HM <- draw(HM)

# **Extract consistent cluster information**
rcl.list <- row_order(HM)  # Ensures extracted clusters match heatmap

# **Save heatmap image**
png(filename = "Sig_Genes_AllGroups_HM_clusters.png", width = 6000, height = 6000, res = 900)
draw(HM)  # Ensures the exact same clustering is saved
dev.off()

# **Convert cluster information into a data frame**
clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(hmtable)[rcl.list[[i]]],
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>% do.call(rbind, .)

# **Map gene symbols to Ensembl IDs**
clu_df$GeneName <- mapIds(org.Hs.eg.db, keys = clu_df$GeneID, keytype = "ENSEMBL", column = "SYMBOL")

# **Save cluster assignments to a text file**
write.table(clu_df, file= "HM_gene_clusters.txt", sep="\t", quote=F, row.names=FALSE)

# Check the number of genes in each cluster from row_order(HM)
cluster_sizes_heatmap <- sapply(rcl.list, length)
print(cluster_sizes_heatmap)

# Check the number of genes in each cluster from the output file
cluster_sizes_file <- table(clu_df$Cluster)
print(cluster_sizes_file)



#venn diagram
venn_list <- list("Distillates and Disposables vs. PG/VG Control" = sig_Dis_vs_Con, 
                  "Juices vs. PG/VG Control" = sig_Ju_vs_Con, 
                  "Distillates and Disposables vs. Juices" = sig_Dis_vs_Ju)

tiff(filename = "Sig_Venn.tiff", width = 2500, height = 2000, res = 300)
ggVennDiagram(venn_list)
dev.off()

# --- 1. Load Dependencies ---
if (!requireNamespace("lobstr", quietly = TRUE)) install.packages("lobstr")
library(gplots)
library(eulerr)
library(lobstr)

# --- 2. Generate Venn and Extract Intersections ---
my_venn <- venn(venn_list)

# Verify structure
print(names(attributes(my_venn)))

# Extract and immediately flatten to base vectors
raw_intersections <- attr(x = my_venn, "intersections")
flat_intersections <- lapply(raw_intersections, as.character) # Forces gene names to character strings

# --- 3. Robust Column Alignment for CSV ---
# This unlist(lapply) combo is the definitive fix for the 'type (list)' error
all_lengths <- unlist(lapply(flat_intersections, length))
max_rows <- max(all_lengths, na.rm = TRUE)

# Build the matrix by padding shorter columns with NA
venn_matrix <- do.call(cbind, lapply(flat_intersections, function(column) {
  # Handle empty intersections safely
  if (length(column) == 0) return(rep(NA, max_rows))
  
  # Pad to max height
  new_col <- column
  length(new_col) <- max_rows
  return(new_col)
}))

# Save to CSV - ready for Google Sheets
write.csv(as.data.frame(venn_matrix), 
          file = "EulerIntersections_VapedD8_3comparisons.csv", 
          row.names = FALSE)

# --- 4. Euler Visualization ---
fit1 <- euler(venn_list)
product_colors <- c("gray40", "white", "lightblue")

png(filename = "Sig_Euler_fit1_3comparisons.png", width = 2000, height = 2000, res = 600)
plot(fit1, 
     legend = list(side = "bottom"),
     quantities = list(type = "counts"),
     fills = product_colors) 
dev.off()

#_______VOLCANOS____________
library(ggrepel)

# Ensure your dataframe has the correct structure
res2$gene <- ifelse(is.na(res2$symbol) | res2$symbol == "", rownames(res2), res2$symbol)

# Define threshold for significance
log2FC_threshold <- 1  # Adjust as needed
pval_threshold <- 0.05

# Add a column for significance
res2$Significance <- "Not Significant"
res2$Significance[res2$logFC > log2FC_threshold & res2$adj.P.Val < pval_threshold] <- "Upregulated"
res2$Significance[res2$logFC < -log2FC_threshold & res2$adj.P.Val < pval_threshold] <- "Downregulated"

# Convert to factor for proper color mapping
res2$Significance <- factor(res2$Significance, levels = c("Upregulated", "Downregulated", "Not Significant"))

# Create the volcano plot
volcano_plot <- ggplot(res2, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
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

res2$padj_num <- as.numeric(res2$adj.P.Val)
top_genes <- res2[order(res2$padj_num), ][1:10, ]  # Adjust number as needed

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
ggsave("JuvsCon_volcano.png", plot = volcano_plot, width = 6.32, height = 4.35, dpi = 900)





#___________________________________________________#
#Making comparisons between each product and the PGVG control
prod <- sapply(snames, function(x) {
  strsplit(x, "_")[[1]][3]  # Split by "_" and get the 2rd part
}) 
print(prod)

prod <- as.factor(prod)
prod

base_names <- gsub("_[0-9]+$", "", snames)
replicate <- factor(base_names)

mm <- model.matrix(~0+prod+batch)

y <- voom(d, mm, plot = T)

#performing duplicateCorrelation to account for tech reps
corfit <- duplicateCorrelation(y, mm, block = replicate)
fit <- lmFit(y, mm, block = replicate, correlation = corfit$consensus.correlation)

contrast.matrix <- makeContrasts(
  A_vs_H = prodA - prodH,
  C_vs_H = prodC - prodH,
  D_vs_H = prodD - prodH,
  E_vs_H = prodE - prodH,
  F_vs_H = prodF - prodH,
  G_vs_H = prodG - prodH,
  levels = colnames(mm)
)

fit_contrast <- contrasts.fit(fit, contrast.matrix)
fit_contrast <- eBayes(fit_contrast)

#Extracting each comparison
top.table_A_vs_H <- topTable(fit_contrast, coef = "A_vs_H", adjust.method = "BH", number = Inf)
top.table_C_vs_H <- topTable(fit_contrast, coef = "C_vs_H", adjust.method = "BH", number = Inf)
top.table_D_vs_H <- topTable(fit_contrast, coef = "D_vs_H", adjust.method = "BH", number = Inf)
top.table_E_vs_H <- topTable(fit_contrast, coef = "E_vs_H", adjust.method = "BH", number = Inf)
top.table_F_vs_H <- topTable(fit_contrast, coef = "F_vs_H", adjust.method = "BH", number = Inf)
top.table_G_vs_H <- topTable(fit_contrast, coef = "G_vs_H", adjust.method = "BH", number = Inf)

# Save or view results
write.csv(top.table_A_vs_H, "A_vs_H.csv")
write.csv(top.table_C_vs_H, "C_vs_H.csv")
write.csv(top.table_D_vs_H, "D_vs_H.csv")
write.csv(top.table_E_vs_H, "E_vs_H.csv")
write.csv(top.table_F_vs_H, "F_vs_H.csv")
write.csv(top.table_G_vs_H, "G_vs_H.csv")

#Making results file of only significant genes for each comparison
#Then making a vector of significant gene probe IDs to use in heatmap
res1 <- top.table_A_vs_H
res2 <- top.table_C_vs_H
res3 <- top.table_D_vs_H
res4 <- top.table_E_vs_H
res5 <- top.table_F_vs_H
res6 <- top.table_G_vs_H

res_sig1 <- res1[which(res1$adj.P.Val < 0.05 & abs(res1$logFC) >= 1), ] 
summary(res_sig1)
sum(res_sig1$logFC > 0, na.rm = TRUE)
sum(res_sig1$logFC < 0, na.rm = TRUE)
res_sig1$symbol <- mapIds(org.Hs.eg.db, keys = res_sig1$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig1, file="A_vs_H_sig.csv")
sig_A_vs_H <- res_sig1$ID  

res_sig2 <- res2[which(res2$adj.P.Val < 0.05 & abs(res2$logFC) >= 1), ] 
summary(res_sig2)
sum(res_sig2$logFC > 0, na.rm = TRUE)
sum(res_sig2$logFC < 0, na.rm = TRUE)
res_sig2$symbol <- mapIds(org.Hs.eg.db, keys = res_sig2$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig2, file="C_vs_H_sig.csv")
sig_C_vs_H <- res_sig2$ID 

res_sig3 <- res3[which(res3$adj.P.Val < 0.05 & abs(res3$logFC) >= 1), ] 
summary(res_sig3)
sum(res_sig3$logFC > 0, na.rm = TRUE)
sum(res_sig3$logFC < 0, na.rm = TRUE)
res_sig3$symbol <- mapIds(org.Hs.eg.db, keys = res_sig3$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig3, file="D_vs_H_sig.csv")
sig_D_vs_H <- res_sig3$ID 

res_sig4 <- res4[which(res3$adj.P.Val < 0.05 & abs(res4$logFC) >= 1), ] 
summary(res_sig4)
sum(res_sig4$logFC > 0, na.rm = TRUE)
sum(res_sig4$logFC < 0, na.rm = TRUE)
res_sig4$symbol <- mapIds(org.Hs.eg.db, keys = res_sig4$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig4, file="E_vs_H_sig.csv")
sig_E_vs_H <- res_sig4$ID 

res_sig5 <- res5[which(res5$adj.P.Val < 0.05 & abs(res5$logFC) >= 1), ] 
summary(res_sig5)
sum(res_sig5$logFC > 0, na.rm = TRUE)
sum(res_sig5$logFC < 0, na.rm = TRUE)
res_sig5$symbol <- mapIds(org.Hs.eg.db, keys = res_sig5$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig5, file="F_vs_H_sig.csv")
sig_F_vs_H <- res_sig5$ID 

res_sig6 <- res6[which(res6$adj.P.Val < 0.05 & abs(res6$logFC) >= 1), ] 
summary(res_sig6)
sum(res_sig6$logFC > 0, na.rm = TRUE)
sum(res_sig6$logFC < 0, na.rm = TRUE)
res_sig6$symbol <- mapIds(org.Hs.eg.db, keys = res_sig6$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res_sig6, file="G_vs_H_sig.csv")
sig_G_vs_H <- res_sig6$ID 

#Filtering for cilia related genes
cilia_genes <- c(
  "Arl4d", "Arl6", "ATMIN", "B9d1", "B9d2", "Bbs1", "Bbs2", "Bbs4", "Bbs7", "Bbs9",
  "C10orf63", "C6orf97", "CAPS", "Cav3", "CCDC", "CCDC103", "Ccdc104", "CCDC114",
  "CCDC151", "CCDC39", "CCDC40", "CCNO", "CEP", "Cep72", "CETN2", "CFAP", "Cfap251",
  "CGI-38", "Cluap1", "collectrin", "CROCC", "DNAAF4", "Dnah", "Dnah11", "Dnah3",
  "DNAH5", "Dnah9", "Dnai1", "DNAI2", "Dnai3", "DNAL1", "Dnali1", "DYDC2", "Dync2h1",
  "Dync2li1", "DYNLL1", "Dynlrb2", "DYX1C1", "Efhc1", "EfhC2", "FKH-2", "FLJ13946",
  "FLJ23577", "FOXJ1", "GHR", "HNF1B", "HSPA1A", "IFT", "Ift172", "Ift20", "Ift46",
  "Ift52", "Ift55", "Ift56", "IFT57", "Ift70b", "Ift80", "IFT81", "Ift88", "Iqca1",
  "IQCE", "Kif9", "LRRC23", "LRRC6", "MCIDAS", "Mks1", "MLF1", "MNS1", "Mucin5A",
  "NEK5", "NME5", "Nphp1", "Nphp4", "NR2F1", "Pkhd1", "RABL4", "RABL5", "Rfx",
  "RIBC2", "ROPN1L", "RSPH", "Rsph4a", "RTDR1", "Spa17", "SPAG6", "Spag6l", "Spata17",
  "SPEF2", "Syne1", "TEKT", "TEKT2", "Tekt4", "Tm4sf1", "Traf3ip1", "Trim28", "Trpv4",
  "Ttc8", "TTLL", "TUBA3", "TUBA4", "TUBB2", "TULP"
)
length(cilia_genes) 

sig_list <- list(res_sig1, res_sig2, res_sig3, res_sig4, res_sig5, res_sig6)
sig_cilia <- lapply(sig_list, function(df) {
  df[grep(paste0("^(", paste(cilia_genes, collapse = "|"), ")"), df$symbol, ignore.case = TRUE), ]
})
names(sig_cilia) <- paste0("sig_cilia", 1:6)
sig_cilia[[1]]
sig_cilia[[2]]
sig_cilia[[3]]
sig_cilia[[4]]
sig_cilia[[5]]
sig_cilia[[6]]

all_sig_genes <- c(sig_A_vs_H, sig_C_vs_H, sig_D_vs_H, sig_E_vs_H, sig_F_vs_H, sig_G_vs_H)
all_sig_genes_unique <- unique(all_sig_genes)

#Making matrix of normalized counts to use in heatmap
normalized_counts <- as.data.frame(y$E)       # Extract normalized counts
rownames(normalized_counts) <- y$genes$ID

base_names <- gsub("_[0-9]+$", "", colnames(normalized_counts))
unique_samples <- unique(base_names)

#average normalized counts for technical replicates
averaged_counts <- matrix(0, nrow = nrow(normalized_counts), ncol = length(unique_samples))
rownames(averaged_counts) <- rownames(normalized_counts)
colnames(averaged_counts) <- unique_samples

for (sample in unique_samples) {
  replicate_cols <- which(base_names == sample)  # Columns corresponding to this sample
  averaged_counts[, sample] <- rowMeans(normalized_counts[, replicate_cols, drop = FALSE])
}
head(averaged_counts)

normcounts <- averaged_counts

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
donors <- c("DD001Rp1", "DD004Op1", "DD008Rp1", "DD021Rp1")
donor_colors <- setNames(RColorBrewer::brewer.pal(4, "Set2"), donors)
coldata <- read.csv('Metadata_notechrep_RawCounts_VapedD8RNASeq.csv', header = TRUE, sep = ",")

library(ComplexHeatmap)

tiff(filename = "Sig_Genes_AllGroups_HM.tiff", width = 2500, height = 2000, res = 300)
Heatmap(hmtable,
        #column_dend_reorder = FALSE,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(
          Product = coldata$Product, 
          Donor = coldata$Donor,
          Type = coldata$Product_type,
          simple_anno_size = unit(.4, "cm"),
          show_annotation_name = FALSE,  
          col = list(Product = c("A" = "#66c2a5", "C" = "#fc8d62", "D" = "#8da0cb", "E" = "#e78ac3", "F" = "#a6d854", "G" = "#ffd92f", "H" = "#e5c494"), 
                     Donor = c("DD001Rp1" = "dodgerblue","DD004Op1" = "forestgreen","DD008Rp1" = "tomato","DD021Rp1" = "gold"),
                     Type = c("Distillate" = "#ffb3b3", "Disposable" = "#c2d9f0", "Juice" = "#b3e0a4", "PGVG" = "#ffcc99")),
          annotation_legend_param = list(Type = list(labels = c("Distillate", "Disposable", "Juice","PGVG"), at = c("Distillate", "Disposable", "Juice","PGVG")))))
dev.off()

#extracting gene lists for specific clusters
HM <- Heatmap(hmtable, km=6)
HM

tiff(filename = "Sig_Genes_AllGroups_HM.tiff", width = 2500, height = 2500, res = 300)
Heatmap(hmtable,
        km = 6,
        #column_dend_reorder = FALSE,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(
          Product = coldata$Product, 
          Donor = coldata$Donor,
          Type = coldata$Product_type,
          simple_anno_size = unit(.4, "cm"),
          show_annotation_name = FALSE,  
          col = list(Product = c("A" = "#66c2a5", "C" = "#fc8d62", "D" = "#8da0cb", "E" = "#e78ac3", "F" = "#a6d854", "G" = "#ffd92f", "H" = "#e5c494"), 
                     Donor = c("DD001Rp1" = "dodgerblue","DD004Op1" = "forestgreen","DD008Rp1" = "tomato","DD021Rp1" = "gold"),
                     Type = c("Distillate" = "#ffb3b3", "Disposable" = "#c2d9f0", "Juice" = "#b3e0a4", "PGVG" = "#ffcc99")),
          annotation_legend_param = list(Type = list(labels = c("Distillate", "Disposable", "Juice","PGVG"), at = c("Distillate", "Disposable", "Juice","PGVG")))))
dev.off()

#VENN
#ggVennDiagram works for up to 7 comparisons, just add more comparisons to the list function below
#using the vectors of significant genes made previously
venn_list <- list("Palmetto Distillate vs. PG/VG Conrol" = sig_A_vs_H, 
                  "NYSW Distillate vs. PG/VG Control" = sig_C_vs_H, 
                  "Palmetto Disposable vs. PG/VG Control" = sig_D_vs_H,
                  "Stiizy Disposable vs. PG/VG Control" = sig_E_vs_H,
                  "ST Strawberry Juice vs. PG/VG Control" = sig_F_vs_H,
                  "ST Green Apple Juice vs. PG/VG Control" = sig_G_vs_H)

tiff(filename = "Sig_Venn.tiff", width = 2500, height = 2000, res = 300)
ggVennDiagram(venn_list)
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
plot(fit1)
fit2 <- euler(venn_list, shape = "ellipse")
fit2
plot(fit2)

product_colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#B3DE69", "#FDB462")

png(filename = "Sig_Euler_fit1.png", width = 4000, height = 4000, res = 600)
plot(fit1, 
     legend = list(side = "bottom"),
     quantities = list(type = c("counts")),
     fills = product_colors)  # Apply colors
dev.off()

tiff(filename = "Sig_Euler_fit2.tiff", width = 2000, height = 2000, res = 300)
plot(fit2, 
     legend = list(side = "bottom"),
     quantities = list(type = c("counts")),
     fills = product_colors)  # Apply colors
dev.off()

#to get gene lists from certain intersections in your Venn diagram
intersect <- intersect(sig_JuicevsApical, sig_OilvsApical) #getting overlapping probe IDs
JuiceOnly <- setdiff(sig_JuicevsApical, sig_OilvsApical)
OilOnly <- setdiff(sig_OilvsApical, sig_JuicevsApical)


intersect <- as.data.frame(intersect) #making it a dataframe
colnames(intersect) = c("GeneID")
intersect$GeneName <- mapIds(org.Hs.eg.db, keys = intersect$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
intersect$GeneName
write.csv(intersect,file = "JuiceandOil_intersection.csv")

JuiceOnly <- as.data.frame(JuiceOnly) #making it a dataframe
colnames(JuiceOnly) = c("GeneID")
JuiceOnly$GeneName <- mapIds(org.Hs.eg.db, keys = JuiceOnly$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
JuiceOnly$GeneName
write.csv(JuiceOnly, file ="JuiceGenesOnly.csv")

OilOnly <- as.data.frame(OilOnly) #making it a dataframe
colnames(OilOnly) = c("GeneID")
OilOnly$GeneName <- mapIds(org.Hs.eg.db, keys = OilOnly$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
OilOnly$GeneName
write.csv(OilOnly,file = "OilGenesOnly.csv")

#_______________________dose response analysis_______________________#
counts <- cts

d0 <- DGEList(counts)# Create DGEList object
d0 <- calcNormFactors(d0)
cutoff <- 4 # genes expressed in at least 10 samples to be included
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

#running linear model for D8THC measured in the vaped products
str(coldata)
head(coldata$D8THC)
rownames(d$counts) <- d$genes$ID

coldata$D8THC <- as.numeric(coldata$D8THC)
design <- model.matrix(~ coldata$D8THC + batch)

y <- voom(d, design, plot = TRUE)

corfit <- duplicateCorrelation(y, design, block = replicate)
fit <- lmFit(y, design, block = replicate, correlation = corfit$consensus.correlation)
fit <- eBayes(fit)
res <- topTable(fit, coef="coldata$D8THC", number=Inf)

res_linD8_sig <- res[which(res$adj.P.Val < 0.05), ] 
sum(res_linD8_sig$logFC > 0, na.rm = TRUE)
sum(res_linD8_sig$logFC < 0, na.rm = TRUE)
res_linD8_sig$symbol <- mapIds(org.Hs.eg.db, keys = res_linD8_sig$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names


median(fit$sigma, na.rm = TRUE)

#visualizing dose response curves for top genes effected by dose
library(reshape2)
top_genes <- rownames(res[res$adj.P.Val < 0.05, ])[1:10]
expr_data <- y$E[rownames(y$E) %in% top_genes, ]
long_expr <- melt(data.frame(D8THC = coldata$D8THC, t(expr_data)), id.vars = "D8THC")

top_genes
gene_symbols <- mapIds(org.Hs.eg.db, keys = top_genes, 
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbols

long_expr$GeneSymbol <- factor(long_expr$variable, levels = names(gene_symbols), labels = gene_symbols)

png(filename = "LinearDoseResponse_Top10Genes_VapedD8_D8.png", width = 2000, height = 1400, res = 300)
ggplot(long_expr, aes(x = D8THC, y = value, color = GeneSymbol)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +  # Linear fit per gene
  theme_minimal() +
  labs(title = "Dose-Response Curves for Top Genes",
       x = "∆8-THC Concentration (nmol/mg)",
       y = "Normalized Expression (log2 CPM)",
       color = "Gene") + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),  # Increase tick label size
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.title = element_text(size = 15)) 
dev.off()



#running linear model for D8THCQ measured in the vaped products
str(coldata)
head(coldata$D8THCQ)
rownames(d$counts) <- d$genes$ID

coldata$D8THCQ <- as.numeric(coldata$D8THCQ)
design <- model.matrix(~ coldata$D8THCQ + batch)

y <- voom(d, design, plot = TRUE)

corfit <- duplicateCorrelation(y, design, block = replicate)
fit <- lmFit(y, design, block = replicate, correlation = corfit$consensus.correlation)
fit <- eBayes(fit)
res <- topTable(fit, coef="coldata$D8THCQ", number=Inf)

res_linD8Q_sig <- res[order(res$adj.P.Val), ][1:500, ]
sum(res_linD8Q_sig$logFC > 0, na.rm = TRUE)
sum(res_linD8Q_sig$logFC < 0, na.rm = TRUE)
res_linD8Q_sig$symbol <- mapIds(org.Hs.eg.db, keys = res_linD8Q_sig$ID, keytype = "ENSEMBL", column = "SYMBOL") #adding gene names

median(fit$sigma, na.rm = TRUE)

#visualizing dose response curves for top genes effected by dose
library(reshape2)
top_genes <- rownames(res[res$adj.P.Val < 0.05, ])[1:10]
expr_data <- y$E[rownames(y$E) %in% top_genes, ]
long_expr <- melt(data.frame(D8THCQ = coldata$D8THCQ, t(expr_data)), id.vars = "D8THCQ")

top_genes
gene_symbols <- mapIds(org.Hs.eg.db, keys = top_genes, 
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbols

long_expr$GeneSymbol <- factor(long_expr$variable, levels = names(gene_symbols), labels = gene_symbols)


png(filename = "LinearDoseResponse_Top10Genes_VapedD8_D8Q.png", width = 2000, height = 1400, res = 300)
ggplot(long_expr, aes(x = D8THCQ, y = value, color = GeneSymbol)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +  # Linear fit per gene
  theme_minimal() +
  labs(title = "Dose-Response Curves for Top Genes",
       x = "∆8-THCQ Concentration (nmol/mg)",
       y = "Normalized Expression (log2 CPM)",
       color = "Gene") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),  # Increase tick label size
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.title = element_text(size = 15)) # Adjust legend title size# Increase axis label size

dev.off()


#making euler to compare d8q reponsive genes to d8q exposure
resd8q_apical <- sig_D8THCQvsVC
resd8q_vaped <- res_linD8Q_sig$ID

venn_list <- list("Genes Altered by ∆8-THCQ" = resd8q_apical, 
                  "Dose Responsive Genes to ∆8-THCQ in Vaped ∆8-THC" = resd8q_vaped) 
                  
fit1 <- euler(venn_list)
fit1
plot(fit1)

png(filename = "d8q_euler.png", width = 2500, height = 2000, res = 600)
plot(fit1, 
     legend = list(side = "bottom"),
     quantities = list(type = c("counts")))
dev.off()

intersect <- intersect(resd8q_apical, resd8q_vaped) #getting overlapping probe IDs
intersect <- as.data.frame(intersect) #making it a dataframe
colnames(intersect) = c("GeneID")
intersect$GeneName <- mapIds(org.Hs.eg.db, keys = intersect$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
intersect$GeneName
write.csv(intersect,file = "d8q17genes.csv")

