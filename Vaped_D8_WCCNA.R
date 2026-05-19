#################################################################################################
#################################################################################################
#### Correlating chemical analysis of ∆8-THC vaped condensates with RNA-seq data from primary human bronchial epithelial cells exposes to vaped ∆8-THC aerosols
#### WCCNA analysis code adapted from Wildfire WGCNA Code (Rager)
#################################################################################################
#################################################################################################

# Clean your working environment
rm(list=ls())

#################################################################################################
#### Setting working directory
#### In RStudio go to Session -> Set Working Directory -> Choose Directory to easily get path
#################################################################################################

#Example
setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/WCCNA Example")


#################################################################################################
#### Installing required R packages
#################################################################################################

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# ---- CRAN packages ----
cran_pkgs <- c(
  "WGCNA",         # network analysis + utilities (numbers2colors, corPvalueStudent)
  "data.table",    # fast data handling
  "tidyverse",     # ggplot2, dplyr, tidyr, readr, etc.
  "reshape2",      # melt/cast (legacy; you use melt)
  "matrixStats",   # colVars used for variance filtering
  "RColorBrewer",  # palettes for plots
  "dynamicTreeCut",# used by WGCNA
  "flashClust",    # fast hierarchical clustering
  "Hmisc",         # general utilities (corr, labeling)
  "circlize",       # colorRamp2 (used by ComplexHeatmap)
  "patchwork"      # combine ggplots with +, /, |
)
to_install_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install_cran)) install.packages(to_install_cran)

# ---- Bioconductor packages ----
bioc_pkgs <- c(
  "DESeq2",         # RNA-seq normalization & DE
  "genefilter",     # filtering helpers
  "AnnotationDbi",  # mapIds framework
  "org.Hs.eg.db",   # human gene annotations (ENSEMBL -> SYMBOL)
  "GO.db",          # GO terms (not strictly required but you listed)
  "preprocessCore", # normalization utils (WGCNA ecosystem)
  "impute",         # imputation (WGCNA ecosystem)
  "ComplexHeatmap"  # heatmaps
)
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) BiocManager::install(to_install_bioc)

#################################################################################################
#### Loading R packages required for this session
#################################################################################################

library(WGCNA)
library(data.table)
library(tidyverse)     # loads ggplot2, dplyr, tidyr, tibble, etc.
library(reshape2)
library(matrixStats)
library(RColorBrewer)
library(dynamicTreeCut)
library(flashClust)
library(Hmisc)

library(DESeq2)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(preprocessCore)
library(impute)

library(ComplexHeatmap)
library(circlize)      
library(grid)
library(patchwork)

# Check that packaes are loaded
sessionInfo()


#################################################################################################
#### Step 1: Loading and organizing chemical data
#################################################################################################

# Read in the data that you want to group, and then eventually correlate to outcomes of interest
# In this case, it's the chemistry data. In genomic-based analyses, this is the -omic profile data
Chemicals = read.csv("Chemistry.csv")
sapply(Chemicals,class) # check classes of each column

# This particular data set includes GC-MS and LC-MS analysis together so the values need to be normalized
# To normalize, we are performing z-score normalization within each chemical row-wise

# Removing chemicals with 0s in all samples
chem_data_filtered <- Chemicals[rowSums(Chemicals[ , -1] != 0) > 0, ]

# Converting data to z-scores
chem_names <- chem_data_filtered[, 1] # Extracting the first column (chemical names) to keep as identifiers
data_matrix <- chem_data_filtered[, -1] # Extracting numeric sample values only for analysis
data_matrix <- as.data.frame(lapply(data_matrix, as.numeric))  # convert all columns to numeric
z_data <- t(scale(t(data_matrix), center = TRUE, scale = TRUE)) # z-score normalize each chemical across samples
z_df <- cbind(Chemical = chem_names, as.data.frame(z_data)) # add chemical names back to normalized data

# Creating a transposed data frame of the chemistry data for future steps
t_Chemicals=as.data.frame(t(z_df[,-c(1)]))
names(t_Chemicals)=z_df$Chemical

# Diagnostic plot
z_melt <- melt(z_df, id.vars = "Chemical", variable.name = "Sample", value.name = "Zscore")
ggplot(z_melt, aes(x = Chemical, y = Zscore)) +
  geom_boxplot(outlier.shape = NA, fill = "darkorange", alpha = 0.6) +
  theme_bw() +
  labs(title = "Z-Score Distribution for Each Chemical",
       x = "Chemical",
       y = "Z-score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
z_df <- cbind(Chemical = chem_names, as.data.frame(z_data))


#################################################################################################
#### Step 2: Loading, filtering, and normalizing transcriptomic data
# Read in the data that represent your outcomes of interest, referred to as "trait" data in WGCNA applications
# Here, they are RNA-seq data, we are starting from raw counts
#################################################################################################

# Loading in data
## Note on sample IDs: A,C,D,E,F,G are different ∆8-THC products and H is the control, the middle number is the shortened donor code, the final number is the technical replicate 
cts = read.csv("RNASeq_Data.csv", header=TRUE, row.names=1) # Loading raw counts
coldata = read.csv("Metadata_RNASeq.csv", header=TRUE, row.names=1) # Loading metadata 
sapply(cts,class) # check classes of each column

# Creating a DESeq2 dataset object from the count matrix and sample metadata (no design formula for now)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 1)

### Collapsing technical replicates
# This strips a trailing _1/_2 
bio_id <- sub("_[0-9]+$", "", colnames(cts)) # map tech reps to their biological sample ID, this drops the _1 or _2

dds_raw   <- DESeqDataSetFromMatrix(countData = cts, # sum raw counts across technical replicates
                                    colData   = coldata,
                                    design    = ~ Product)

dds_coll  <- collapseReplicates(dds_raw, groupby = bio_id, run = colnames(cts))
cts_coll  <- counts(dds_coll)
coldata_coll <- as.data.frame(colData(dds_coll))

#### Descriptive statistics to determine how to filter genes
#Around 2k genes is a good number for WCCNA

med_ct <- apply(cts_coll, 1, median, na.rm = TRUE) # per-gene median count after collapsing tech reps
quantile(med_ct, probs = c(0.25, 0.5, 0.75)) # quartiles

# Plot histogram of log10-transformed median counts and mark cutoff lines for 25 (blue) and 50 (red)
hist(log10(med_ct + 1), breaks = 100,
     main = "Log10 distribution of median counts per gene",
     xlab = "log10(median raw count + 1)")
abline(v = log10(25), col = "blue", lwd = 2, lty = 2) # Genes with median raw count of 25
abline(v = log10(50), col = "red", lwd = 2, lty = 2) # Genes with median raw count of 50
#many genes have low or no expression

# Number of genes kept at each threshold
sum(med_ct >= 25)  # genes kept at 25
sum(med_ct >= 50)  # genes kept at 50

cts_f <- counts(dds_coll)[med_ct >= 25, , drop = FALSE] #using median count >=25 (16790 genes)

#### Next going to cut down genes even more
# Filtering by genes altered in more than one exposure vs. the control (H)

coldata_coll$Product <- factor(coldata_coll$Product)        # must exist in your metadata
coldata_coll$Product <- relevel(coldata_coll$Product, "H")  # set control (H)
stopifnot(all(colnames(cts_f) == rownames(coldata_coll)))   # ensure sample order matches
dds_de  <- DESeqDataSetFromMatrix(countData = cts_f, colData = coldata_coll, design = ~ Product) # create DESeq2 dataset using filtered counts and exposure groups for contrasts
dds_de  <- DESeq(dds_de) # run DESeq2

# run contrasts for each exposure vs control
alph <- 0.05
lfc_cut <- log2(1.5)   # ~0.585 = 1.5x fold change cutoff
expos   <- setdiff(levels(colData(dds_de)$Product), "H")

res_list <- lapply(expos, function(x) {
  res <- results(dds_de, contrast = c("Product", x, "H"))
  sig <- (res$padj < alph) & (abs(res$log2FoldChange) >= lfc_cut)
  sig[is.na(sig)] <- FALSE
  sig
})
sig_mat <- do.call(cbind, res_list) # logical matrix: rows = genes, columns = exposures; TRUE if gene is significant in that exposure vs H
colnames(sig_mat) <- expos

sig_count <- rowSums(sig_mat) # for each gene, how many exposures it’s significant in
table(sig_count)
sum(sig_count >= 3)   # adjust k as needed, this counts genes are significant in at least 3 exposures, gets us to 2073

rownames(sig_mat) <- rownames(cts_f)   # attach gene IDs
sig_count <- rowSums(sig_mat)
robust_genes <- names(sig_count)[sig_count >= 3]  # genes significant in ≥3 exposures

#### Normalizing
vsd <- vst(dds_de, blind = TRUE) # Variance stabilizing transformation
expr_vst <- assay(vsd)[robust_genes, , drop = FALSE]   # genes x samples


#################################################################################################
#### Step 3: Setting up WCCNA and initial visualization of data
## Chemistry data will be tied to donor and exposure combination
## Specifically making a dendogram of how trait (toxicity) data relate to sample chemistry data
#################################################################################################

Traits <- as.data.frame(t(expr_vst)) # Putting normalized, filtered RNA-seq seq data into Traits

# Create product label (A/C/D/E/F/G) from sample IDs like "A_21"
prod_from_traits <- sub("_.*", "", rownames(Traits))                 

# Expand product-level chemistry to donor-level by repeating each product’s row
chem_expanded <- t_Chemicals[prod_from_traits, , drop = FALSE]    
rownames(chem_expanded) <- rownames(Traits)
stopifnot(identical(rownames(Traits), rownames(chem_expanded)))          # alignment check

# Remove control H from both matrices (no chemistry measured for H)
keep <- prod_from_traits != "H"
Traits        <- Traits[keep, , drop = FALSE]
chem_expanded <- chem_expanded[keep, , drop = FALSE]
stopifnot(identical(rownames(Traits), rownames(chem_expanded)))

# Cluster samples based on chemistry
sampleTree <- hclust(dist(chem_expanded), method = "average")

# Align RNA-seq traits to chemistry rows (same samples, same order)
Traits_matched <- Traits[rownames(chem_expanded), , drop = FALSE]        # align to chem_expanded

# select top 100 most variable genes across samples
gene_vars <- colVars(as.matrix(Traits_matched)) # variance per gene (columns)
top_genes <- order(gene_vars, decreasing = TRUE)[1:100] # indices of top-variance genes
top_Traits <- Traits_matched[, top_genes] # subset to top 100 genes

# hierarchical clustering of samples (the donor × product combinations), based on their expanded chemistry profiles.
traitColors <- numbers2colors(scale(top_Traits), signed = FALSE)
plotDendroAndColors(
  dendro = sampleTree,
  colors = traitColors,
  groupLabels = colnames(top_Traits),
  addGuide = FALSE,
  guideAll = FALSE, 
  guideCount = 50,
  guideHang = 0.5, 
  cex.colorLabels = 0.6,
  cex.dendroLabels = 0.6, 
  cex.rowText = 0.6,
  main = "Sample dendrogram and top 100 RNA-seq genes"
)

# Exporting organized data needed for next steps (chemistry + RNA-seq traits aligned by donor-product sample)
save(chem_expanded, Traits_matched, file = "Chemistry_RNAseq_dataInput.RData")


#################################################################################################
#### Step 4: Chemistry network construction and module detection (WCCNA)
## Build a weighted network of chemicals, identify co-modulated chemical modules,
## and export module eigengenes + assignments for downstream association with RNA-seq traits
#################################################################################################

# Load the aligned chemistry (samples x chemicals) and RNA-seq traits saved earlier
lnames = load(file = "Chemistry_RNAseq_dataInput.RData")
lnames
nrow(chem_expanded); ncol(chem_expanded)

# Choose soft-thresholding power for the chemistry network
# pickSoftThreshold evaluates how well different powers approximate a scale-free topology
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(chem_expanded, powerVector = powers, verbose = 5)
#With only ~25 chemicals, the fit (SFT.R.sq) is pretty poor at all powers
#picking a low/moderate power (e.g., 6–9) where mean connectivity isn’t too high but modules can still be separated. 9 seems reasonable.

# Plot scale-free fit (left) and mean connectivity (right) across candidate powers
par(mfrow = c(1,2)); cex1 <- 0.9
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit (signed R^2)",
     type="n", main="Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h = 0.90, col = "red")  # target R^2 ~0.9

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean connectivity",
     type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, cex=cex1, col="red")

# Use Pearson correlation inside WGCNA (wrapper keeps defaults explicit)
cor <- function(x, y = NULL, use = "everything",
                method = "pearson", weights.x = NULL, weights.y = NULL, cosine = FALSE) {
  stats::cor(x, y, use = use, method = method)
}

# Build the chemistry co-expression network and detect chemical modules
# networkType: unsigned behavior (negative and positive correlations of same magnitude both count)
# TOMType: keep topology consistent with the unsigned adjacency
# minModuleSize: allow small chemical modules (we only have ~25 chemicals)
# mergeCutHeight: low value keeps modules distinct; raise to merge similar modules
net = blockwiseModules(chem_expanded,
                       power = 9,
                       TOMType = "unsigned",
                       minModuleSize = 3,
                       mergeCutHeight = 0.1,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ChemicalsTOM",
                       verbose = 3)
# Note: if you wanted only positively co-modulated chemicals, use networkType="signed", TOMType="signed".


# Convert numeric module labels to colors for visualization
mergedColors = labels2colors(net$colors)

# Plot chemical dendrogram with module color bars underneath 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Store key outputs for downstream association with RNA-seq traits
moduleLabels = net$colors            # per-chemical module assignment (numeric)
moduleColors = labels2colors(net$colors)  # same assignments as colors
MEs = net$MEs                        # module eigengenes (one value per sample per module)
geneTree = net$dendrograms[[1]]      # chemistry dendrogram (naming kept for compatibility)

# Save the key network construction results (module eigengenes, module assignments, dendrogram)
# into an .RData file so you can reload later without rebuilding the whole network
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Chemicals-networkConstruction-auto.RData")

# Export module eigengenes (one representative value per module per sample)
# Useful for downstream correlation with traits (e.g., RNA-seq data) or external analysis
write.csv(MEs, file = "MEs_AggregateValues.csv")

# Export a lookup table of chemical-to-module assignments
# moduleLabels = numeric module ID
# moduleColors = human-friendly color label for plotting
moduleLabels_df <- as.data.frame(moduleLabels)
moduleColors_df <- as.data.frame(moduleColors)
GeneModuleAssignments <- cbind(moduleLabels_df, moduleColors_df)
write.csv(GeneModuleAssignments, file = "GeneModuleAssignments.csv")


#################################################################################################
####Step 5: Relating modules to trait information
## Goal: correlate each chemical module eigengene (rows) with each RNA-seq gene (columns),
## then making plots to explore the data
#################################################################################################

# Load the aligned chemistry × sample data and RNA-seq traits from Step 3
lnames = load(file = "Chemistry_RNAseq_dataInput.RData");

# Load the module assignments and eigengenes from Step 4 (network construction)
lnames = load(file = "Chemicals-networkConstruction-auto.RData");
lnames

# Define dimensions of the chemistry data (needed for correlation tests)
nChemicals = ncol(chem_expanded);
nSamples = nrow(chem_expanded);

# Recalculate module eigengenes (MEs = first principal component of each module)
# These are the single "summary profiles" for each chemical module across samples.
MEs0 = moduleEigengenes(chem_expanded, moduleColors)$eigengenes
MEs = orderMEs(MEs0) # reorder columns so similar modules are adjacent

# Correlate module eigengenes (rows) with RNA-seq traits (columns = genes kept for WCCNA)
# moduleTraitCor   = correlation coefficients (Pearson r)
# moduleTraitPvalue = p-values for each correlation (nSamples is used for the test)
moduleTraitCor = cor(MEs, Traits_matched, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Format correlation and p-value together so both appear in each heatmap cell, e.g. "0.65 (0.01)"
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3)); # widen plot margins so labels fit

# ---- Figure 1: Correlation Matrix (limiting to 30 genes) ----
# Step 1: Compute absolute correlations
absCor <- abs(cor(MEs, Traits_matched, use = "p"))

# Step 2: keep the top N genes whose max correlation (with any module) is largest
topN <- min(30, ncol(Traits_matched))
topTraitIndices <- order(apply(absCor, 2, max, na.rm = TRUE), decreasing = TRUE)[1:topN]
Traits_top <- Traits_matched[, topTraitIndices, drop = FALSE]

# Step 3: recompute correlations and p-values on this reduced set (for the figure)
moduleTraitCor <- cor(MEs, Traits_top, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Step 4: re-make the cell text "r (p)"
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Step 5: make the x-axis readable by mapping Ensembl IDs → gene symbols
library(org.Hs.eg.db)
ensembl_ids <- colnames(Traits_top)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
# Use the gene symbol when available; otherwise fall back to the Ensembl ID
xLabels_named <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)

# Step 6: Plot
png("ModuleTraitHeatmap_topTraits.png", width = 2400, height = 1600, res = 300)

# Make and save the heatmap figure
# This shows:
#   * rows = module eigengenes (one row per chemical module)
#   * columns = top RNA-seq genes (most strongly associated with any module)
#   * colors = correlation strength and direction (blue = negative, red = positive)
#   * numbers = "r" on first line; "(p)" on second line for significance
par(mar = c(12, 14, 4, 4))  # Adjust margin for better axis text
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = xLabels_named,  # << use renamed gene labels
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 0.3,
               zlim = c(-1, 1),
               main = "Module–Trait Relationships (Top Traits)")

dev.off()
   

# ---- Figure 2: Heatmap of chemical abundances (Z-scores) across samples, grouped by module and product ----
# Rows = chemicals, grouped (row_split) by their WGCNA module assignment and annotated by module color
# Columns = samples (donor × product), split by product (A, C, D, E, F, G)
# Fill = Z-scored abundance values (blue = low, red = high)

# Build the chemicals x samples matrix from expanded chemistry
z_matrix <- t(as.matrix(chem_expanded))  # rows = chemicals, cols = samples (e.g., "A_21")

# Make sure moduleColors is named by chemical and aligns to rows of z_matrix
names(moduleColors) <- colnames(chem_expanded)
z_matrix <- z_matrix[names(moduleColors), , drop = FALSE]
row_modules <- moduleColors[rownames(z_matrix)]

# Column split by product (A/C/D/E/F/G)
col_split <- factor(sub("_.*", "", colnames(z_matrix)), levels = c("A","C","D","E","F","G"))

# Row annotation by module color
row_anno <- ComplexHeatmap::rowAnnotation(
  Module = row_modules,
  col = list(Module = structure(unique(row_modules), names = unique(row_modules))),
  show_legend = TRUE
)

# Build heatmap object
ht <- ComplexHeatmap::Heatmap(
  z_matrix,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_split = factor(row_modules, levels = unique(row_modules)),
  left_annotation = row_anno,
  column_split = col_split,
  column_names_rot = 90,
  column_title = "Samples (grouped by product)",
  row_title = "Chemicals (grouped by module)"
)

png("Chemicals_by_Sample_Heatmap.png", width = 2400, height = 1600, res = 300)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# ---- Figure 3: Line plots of chemical abundances by WGCNA module ---- 
# Each facet = one WGCNA module (colored panel)
# X-axis = samples (donor × product)
# Y-axis = Z-scored chemical abundance
# Each line = one chemical belonging to that module
#
# Purpose: This figure shows the within-sample variability of chemicals in each module and how
#          abundance patterns differ across samples and products. Chemicals in the same module 
#          are expected to track together, highlighting module coherence.

# 1. Create lookup of each chemical → its assigned module color
module_df <- data.frame(
  gene_id = colnames(chem_expanded),   # chemical names
  colors = moduleColors,               # module colors from WGCNA
  stringsAsFactors = FALSE
)

# 2. Build long-format dataframe of chemical abundances with module assignments
submod_df <- chem_expanded %>%                            # chemistry matrix: samples × chemicals
  rownames_to_column("Sample") %>%                        # keep sample IDs for plotting
  pivot_longer(-Sample, names_to = "gene_id", values_to = "Zscore") %>%  # long format
  left_join(module_df, by = "gene_id") %>%                # add module assignment (colors)
  filter(!is.na(colors))                                  # keep only chemicals assigned to modules

# 3. Plot: one line per chemical, grouped by module color
ggplot(submod_df, aes(x = Sample, y = Zscore, group = gene_id)) +
  geom_line(aes(color = colors), alpha = 0.3) +           # faint lines = individual chemicals
  facet_wrap(~ colors, scales = "free_y") +               # separate panels by module
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    x = "Sample (donor × product)",
    y = "Z-scored chemical abundance",
    title = "Chemical profiles by module"
  )

# ---- Figure 4: Summary Figure - Product-level chemistry + Top kME genes (like Kim 2022 Fig. 1) ----
## Left: chemicals × products (6 cols; product names via sample_map)
## Right: chemicals × genes (top 5 “item-vector” genes per module)
#For each module, the top 5 genes were selected based on the highest absolute kME values (module membership), so that both positively and negatively associated hub genes were considered

# Pretty names for products (columns on left panel)
sample_map <- c(
  A = "Palmetto Distillate",
  C = "NYSW Distillate",
  D = "Palmetto Disposable",
  E = "Stiiizy Disposable",
  F = "ST Strawberry Juice",
  G = "ST Green Apple Juice"
)

# 1) kME (item vectors) and top genes per module
ME_cols     <- grep("^ME", colnames(MEs), value = TRUE)
MEs_numeric <- MEs[, ME_cols, drop = FALSE]
kME         <- cor(MEs_numeric, Traits_matched, use = "p")        # rows=modules (ME*), cols=genes
rownames(kME) <- sub("^ME", "", rownames(kME))                    # "MEblue" -> "blue"
modules     <- rownames(kME)

# Top 5 genes by |kME| in each module
top5_genes <- lapply(modules, function(mod) {
  v <- kME[mod, ]
  names(sort(abs(v), decreasing = TRUE))[1:5]
})
names(top5_genes) <- modules
selected_genes <- unique(unlist(top5_genes))

# Map Ensembl -> SYMBOL for nicer labels
sym_map <- mapIds(org.Hs.eg.db, keys = selected_genes,
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_labels <- ifelse(is.na(sym_map) | sym_map == "", selected_genes, sym_map)

# 2) Correlate chemicals with selected genes (right panel) 
# chem_expanded: 24×K (samples × chems). We need chemicals × genes.
cor_mat <- cor(chem_expanded, Traits_matched[, selected_genes, drop = FALSE], use = "p")
# cor() above returns K×G already (rows = chemicals, cols = genes)
colnames(cor_mat) <- gene_labels

# 3) Product-level chemistry (left panel; 6 columns)
# Collapse the 24 sample rows in chem_expanded to 6 products by averaging donors
chem_prod <- chem_expanded %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::mutate(product = sub("_.*", "", sample)) %>%
  dplyr::select(-sample) %>%
  dplyr::group_by(product) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop")  # use median if you prefer

# Make it a products × chemicals matrix with product codes as rownames (A,C,D,E,F,G)
chem_prod <- as.data.frame(chem_prod)
rownames(chem_prod) <- chem_prod$product
chem_prod$product <- NULL

# reorder products and apply names for columns in the heatmap
prod_order <- c("A","C","D","E","F","G")
chem_prod  <- chem_prod[prod_order, , drop = FALSE]

# Transpose to chemicals × products for the left heatmap
z_mat_prod <- t(as.matrix(chem_prod))                  # rows = chemicals, cols = products
colnames(z_mat_prod) <- sample_map[colnames(z_mat_prod)]

# 4) Align row (chemical) order across both panels by module
# Make sure moduleColors is named by chemical
names(moduleColors) <- colnames(chem_expanded)                     # chemicals vector

# Compute a consistent chemical order (by module color, then name)
chem_order <- order(moduleColors[rownames(z_mat_prod)], rownames(z_mat_prod))

# Reorder both matrices and the module color vector
z_mat_prod  <- z_mat_prod[chem_order, , drop = FALSE]
cor_mat     <- cor_mat[rownames(z_mat_prod), , drop = FALSE]
row_modules <- moduleColors[rownames(z_mat_prod)]

# 5) Annotations
# Row module color bar (shared across panels)
row_anno <- rowAnnotation(
  Module = factor(row_modules, levels = unique(row_modules)),
  col    = list(Module = setNames(unique(row_modules), unique(row_modules))),
  show_legend = TRUE
)

# Optional chemical IDs strip (left side)
chem_ids <- seq_len(nrow(z_mat_prod))
id_anno <- rowAnnotation(
  ID = anno_text(chem_ids, gp = gpar(fontsize = 8)),
  width = unit(4, "mm"),
  show_legend = FALSE
)

# Top annotation on RIGHT: gene’s module (the module it loads highest on)
gene_to_module <- sapply(colnames(cor_mat), function(label) {
  ens <- names(gene_labels)[gene_labels == label]
  hit <- names(Filter(function(v) ens %in% v, top5_genes))
  if (length(hit)) hit[1] else NA
})
col_anno <- HeatmapAnnotation(
  Module = factor(gene_to_module, levels = modules),
  col    = list(Module = setNames(unique(row_modules), unique(row_modules))),
  show_legend = TRUE
)
row_anno <- rowAnnotation(
  Module = factor(row_modules, levels = unique(row_modules)),
  col = list(Module = setNames(unique(row_modules), unique(row_modules))),
  show_legend = FALSE   # hide this legend
)

# 6) Build heatmaps
# LEFT: chemicals × products
ht_left <- Heatmap(
  z_mat_prod,
  name = "Abundance\n(Z-score)",
  col = colorRamp2(c(min(z_mat_prod, na.rm = TRUE), max(z_mat_prod, na.rm = TRUE)),
                   c("white", "darkgreen")),  # single-hue scale
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = moduleColors[chem_order],   # same split you were using
  left_annotation  = row_anno,         # keep your existing annotations
  right_annotation = id_anno,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,                  # <<< vertical labels
  column_names_gp  = gpar(fontsize = 9),  # smaller text
  column_title = "Products",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title = "Chemical",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# RIGHT: chemicals × selected genes 
ht_right <- Heatmap(
  cor_mat,
  name = "Correlation",
  col = colorRamp2(c(min(cor_mat, na.rm = TRUE), 0, max(cor_mat, na.rm = TRUE)),
                   c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = moduleColors[chem_order],
  top_annotation = col_anno,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,                  # <<< vertical gene labels
  column_names_gp  = gpar(fontsize = 9),  # smaller text
  column_title = "Top 5 Genes per Module",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# 7) Draw & save
png("SplitHeatmap_ProductChem_and_TopKMEGenes.png", width = 2800, height = 1700, res = 320)
draw(ht_left + ht_right,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 15), "mm"))
dev.off()

# Get chemical names corresponding to the numbers in the heatmap
chem_id_lookup <- data.frame(
  ID       = seq_len(nrow(z_mat_prod)),
  Chemical = rownames(z_mat_prod),
  Module   = row_modules,                  # moduleColors aligned to z_mat_prod
  Abbrev   = abbreviate(rownames(z_mat_prod), minlength = 6, strict = TRUE),
  stringsAsFactors = FALSE
)

# Look and/or save
head(chem_id_lookup, 10)
write.csv(chem_id_lookup, "Heatmap_ChemID_lookup.csv", row.names = FALSE)


# ---- Figure 5: Top-5 genes (kME) and module chemicals over samples ----

# --- Helpers ---
.abbrev15 <- function(x) ifelse(nchar(x) > 15, paste0(substr(x, 1, 15), "..."), x)

# Map eigengene line color to module name
eig_color_for <- function(mod) {
  switch(tolower(mod),
         "blue"      = "blue",
         "brown"     = "brown",
         "turquoise" = "turquoise",
         "yellow"    = "yellow",
         "grey"      = "grey",
         "black")
}

# Function: stacked plot (top genes; bottom chemicals + eigengene) for a module 
make_module_stack_plot <- function(mod, offset_map = NULL) {
  mod <- tolower(mod)
  me_col <- paste0("ME", mod)
  
  # 1) Order samples by the (sign-locked) module eigengene
  stopifnot(me_col %in% colnames(MEs))
  me_vec <- as.numeric(MEs[, me_col])
  
  stopifnot(mod %in% rownames(kME))
  top_pos_gene <- names(sort(kME[mod, ], decreasing = TRUE))[1]
  if (cor(me_vec, Traits_matched[, top_pos_gene], use = "p") < 0) me_vec <- -me_vec
  
  ord <- order(me_vec)
  sample_order <- rownames(MEs)[ord]
  
  # 2) TOP: top 5 genes by |kME|
  genes_ens <- names(sort(abs(kME[mod, ]), decreasing = TRUE))[1:5]
  
  sym_map <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes_ens,
                                   column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  gene_labels <- ifelse(is.na(sym_map) | sym_map == "", genes_ens, sym_map)
  names(gene_labels) <- genes_ens
  
  expr_top <- scale(Traits_matched[, genes_ens, drop = FALSE]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    mutate(sample = factor(sample, levels = sample_order)) %>%
    arrange(sample)
  
  colnames(expr_top)[match(genes_ens, colnames(expr_top))] <- gene_labels[genes_ens]
  
  expr_long <- expr_top %>%
    tidyr::pivot_longer(cols = -sample, names_to = "Gene", values_to = "Z")
  
  gene_plot <- ggplot(expr_long, aes(sample, Z, color = Gene, group = Gene)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = setNames(brewer.pal(5, "Set1"), unique(expr_long$Gene))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(title = paste0("Module ", mod, ": top 5 genes (kME)"),
         x = NULL, y = "Gene expression (Z-score)") +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.key.height = unit(8, "pt")
    )
  
  # 3) BOTTOM: all chemicals in the module + eigengene line
  chem_in_mod <- names(moduleColors)[tolower(moduleColors) == mod]
  if (length(chem_in_mod) == 0) {
    warning(sprintf("No chemicals found for module '%s'.", mod))
    return(NULL)
  }
  
  chem_long <- chem_expanded[, chem_in_mod, drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    mutate(sample = factor(sample, levels = sample_order)) %>%
    arrange(sample) %>%
    pivot_longer(cols = -sample, names_to = "Chemical", values_to = "Z")
  
  # 15-char legend labels
  chem_long$ChemLab <- .abbrev15(chem_long$Chemical)
  
  # Optional small offsets for specific overlapping series
  if (!is.null(offset_map)) {
    chem_long <- chem_long %>%
      mutate(Z_offset = ifelse(Chemical %in% names(offset_map),
                               Z + offset_map[Chemical], Z))
  } else {
    chem_long$Z_offset <- chem_long$Z
  }
  
  # Eigengene (scaled) overlay
  eig_label <- paste0(mod, " module eigengene")
  eig_df <- data.frame(
    sample   = factor(sample_order, levels = sample_order),
    Z_offset = as.numeric(scale(me_vec[ord])),
    ChemLab  = eig_label
  )
  
  chem_levels <- sort(unique(chem_long$ChemLab))
  all_levels  <- c(eig_label, chem_levels)
  
  pal <- RColorBrewer::brewer.pal(max(3, length(chem_levels)), "Set2")[seq_along(chem_levels)]
  names(pal) <- chem_levels
  col_vec <- c(setNames(eig_color_for(mod), eig_label), pal)
  
  chem_plot <- ggplot(chem_long, aes(sample, Z_offset, color = ChemLab, group = Chemical)) +
    geom_line(linewidth = 0.9, alpha = 0.9) +
    geom_line(data = eig_df,
              aes(x = sample, y = Z_offset, color = ChemLab, group = 1),
              linewidth = 1.8, alpha = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = col_vec, breaks = all_levels, drop = FALSE) +
    labs(title = paste0("Module ", mod, ": member chemicals + eigengene"),
         x = "Samples (ordered by module eigengene)",
         y = "Chemical abundance (Z-score)", color = NULL) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.key.height = unit(8, "pt"),
      legend.text = element_text(size = 9)
    ) +
    guides(color = guide_legend(ncol = 1))
  
  # Stack
  gene_plot / chem_plot + plot_layout(heights = c(0.9, 1.3))
}

# Build the five modules (blue, brown, turquoise, yellow, grey)
mods <- c("blue", "brown", "turquoise", "yellow", "grey")

# (Optional) tiny offsets for known overlaps; leave empty by default
offsets <- c("2,3,4-Trimethylhexane" = 0.05)  # example used for blue; ignored for others

plots <- lapply(mods, function(m) {
  if (m == "blue") {
    make_module_stack_plot(m, offset_map = offsets)
  } else {
    make_module_stack_plot(m)
  }
})

# Drop any NULLs (in case a module has no chemicals)
plots <- Filter(Negate(is.null), plots)

# Combine into a multi-panel (2 rows x 3 cols works nicely for 5 plots)
multi_5 <- wrap_plots(plots, ncol = 3)

# Also save individual module plots
out_dir <- "module_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (m in mods) {
  p <- if (m == "blue") {
    make_module_stack_plot(m, offset_map = offsets)
  } else {
    make_module_stack_plot(m)
  }
  
  if (!is.null(p)) {
    ggsave(file.path(out_dir, sprintf("Module_%s_genes_and_chemicals.png", m)),
           p, width = 12, height = 8, dpi = 300)
  } else {
    message(sprintf("Skipping '%s' (no chemicals found).", m))
  }
}

# Save
ggsave("Modules_blue_brown_turquoise_yellow_grey.png", multi_5,
       width = 21, height = 13, dpi = 300)



#################################################################################################
#### New: Module-specific split heatmaps (Fig 4 variation)
#### - One figure per module (blue, brown, turquoise, yellow; excludes grey)
#### - Left panel: chemicals × products (Abundance, green scale), product names on bottom
#### - Far-left: numeric compound IDs
#### - Right panel: chemicals × top 10 genes by |kME| for that module (Correlation, blue-white-red)
#### - Top titles: "Compound abundance" (left), "Gene correlation with compound" (right)
#### - Final: tile four module figures into a single 2×2 PNG
#################################################################################################

# If not already installed:
if (!"gridExtra" %in% rownames(installed.packages())) install.packages("gridExtra")

library(gridExtra)
library(grid)

# Safety checks for required objects produced earlier
stopifnot(exists("chem_expanded"), exists("moduleColors"), exists("MEs"),
          exists("Traits_matched"), exists("kME"), exists("sample_map"))

# Rebuild product-level chemistry matrix (chemicals × products), as in Fig 4
chem_prod <- chem_expanded %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::mutate(product = sub("_.*", "", sample)) %>%
  dplyr::select(-sample) %>%
  dplyr::group_by(product) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
  as.data.frame()
rownames(chem_prod) <- chem_prod$product
chem_prod$product <- NULL
prod_order <- c("A","C","D","E","F","G")
chem_prod  <- chem_prod[prod_order, , drop = FALSE]
z_mat_prod_all <- t(as.matrix(chem_prod))                          # rows = chemicals, cols = products
colnames(z_mat_prod_all) <- sample_map[colnames(z_mat_prod_all)]   # pretty names

# Helper to build one module’s combined heatmap (labels flush as row names, no left annotation)
make_module_htlist <- function(mod_color) {
  chems_mod <- names(moduleColors)[tolower(moduleColors) == tolower(mod_color)]
  if (!length(chems_mod)) return(NULL)
  
  # Left: abundance by product (subset for this module)
  z_mat_prod <- z_mat_prod_all[chems_mod, , drop = FALSE]
  
  # Truncated row labels (first 20 chars + …), printed as row names (flush to matrix)
  chem_labels <- ifelse(nchar(rownames(z_mat_prod)) > 20,
                        paste0(substr(rownames(z_mat_prod), 1, 20), "…"),
                        rownames(z_mat_prod))
  
  # Right: correlations with top 10 |kME| genes
  this_mod <- tolower(mod_color)
  stopifnot(this_mod %in% rownames(kME))
  top10_genes <- names(sort(abs(kME[this_mod, ]), decreasing = TRUE))[1:10]
  sym_map10 <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = top10_genes,
                                     column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  gene_labels10 <- ifelse(is.na(sym_map10) | sym_map10 == "", top10_genes, sym_map10)
  names(gene_labels10) <- top10_genes
  
  cor_mat <- cor(chem_expanded[, chems_mod, drop = FALSE],
                 Traits_matched[, top10_genes, drop = FALSE], use = "p")
  colnames(cor_mat) <- gene_labels10
  cor_mat <- cor_mat[rownames(z_mat_prod), , drop = FALSE]
  
  # --- Square cell sizing ---
  cell_mm <- 4  # choose your square size in mm (try 3–6)
  left_w  <- unit(ncol(z_mat_prod) * cell_mm, "mm")
  left_h  <- unit(nrow(z_mat_prod) * cell_mm, "mm")
  
  right_w <- unit(ncol(cor_mat)   * cell_mm, "mm")
  right_h <- unit(nrow(cor_mat)   * cell_mm, "mm")
  
  # Color scales
  green_scale <- circlize::colorRamp2(
    c(min(z_mat_prod, na.rm = TRUE), max(z_mat_prod, na.rm = TRUE)),
    c("white", "darkgreen")
  )
  bwrd <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Left heatmap: use row labels (no left_annotation) so labels sit tight against tiles
  ht_left <- Heatmap(
    z_mat_prod,
    name = "Abundance\n(Z-score)",
    col = green_scale,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE,
    row_labels = chem_labels,
    row_names_gp = gpar(fontsize = 9),
    row_names_max_width = max_text_width(chem_labels, gp = gpar(fontsize = 8)),
    row_names_side = "left",
    show_column_names = TRUE, column_names_rot = 90,
    column_title = "Compound Abundance",
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                labels_gp = gpar(fontsize = 8)),
    width  = left_w,   # << makes cells square
    height = left_h,   # <<
    border = TRUE
  )
  
  # Right heatmap
  ht_right <- Heatmap(
    cor_mat,
    name = "Correlation",
    col = bwrd,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE, column_names_rot = 90,
    column_title = "Gene Correlation with Compounds",
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                labels_gp = gpar(fontsize = 8)),
    width  = right_w,  # << makes cells square
    height = right_h,  # <<
    border = TRUE
  )
  
  # Return combined HeatmapList; gap small so panels sit close
  (ht_left + ht_right)
}

# Build per-module HeatmapLists (no grabbing)
mods4 <- c("blue", "brown", "turquoise", "yellow")
htlists <- lapply(mods4, make_module_htlist)
names(htlists) <- mods4
htlists <- Filter(Negate(is.null), htlists)

# --- Save individual PNGs (draw directly to the device) ---
out_dir_modsplit <- "module_split_heatmaps"
if (!dir.exists(out_dir_modsplit)) dir.create(out_dir_modsplit, recursive = TRUE)

for (m in names(htlists)) {
  png(file.path(out_dir_modsplit, sprintf("Module_%s_split_heatmap.png", m)),
      width = 1800, height = 1200, res = 300)
  grid.newpage()
  ComplexHeatmap::draw(
    htlists[[m]],
    ht_gap = unit(2, "mm"),
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(5, 8, 5, 8), "mm"),
    newpage = FALSE
  )
  grid.text(paste0("Module: ", tools::toTitleCase(m)),
            x = unit(0.02, "npc"), y = unit(0.98, "npc"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold"))
  dev.off()
}

# --- Save tiled 2×2 composite on one page (use viewports; no grob capture) ---
png(file.path(out_dir_modsplit, "Modules_blue_brown_turquoise_yellow_split_heatmaps.png"),
    width = 3200, height = 2200, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))

pos <- matrix(c(1,1, 1,2, 2,1, 2,2), byrow = TRUE, ncol = 2)
mods_in_order <- c("blue","brown","turquoise","yellow")

for (i in seq_along(mods_in_order)) {
  m <- mods_in_order[i]
  if (!m %in% names(htlists)) next
  r <- pos[i,1]; c <- pos[i,2]
  vp <- viewport(layout.pos.row = r, layout.pos.col = c)
  pushViewport(vp)
  ComplexHeatmap::draw(
    htlists[[m]],
    ht_gap = unit(2, "mm"),
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(5, 8, 5, 8), "mm"),
    newpage = FALSE
  )
  grid.text(paste0("Module: ", tools::toTitleCase(m)),
            x = unit(0.02, "npc"), y = unit(0.98, "npc"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold"))
  upViewport()
}
dev.off()

# --- Optional: per-module compound lookup files (unchanged logic) ---
for (m in mods4) {
  chems_mod <- names(moduleColors)[tolower(moduleColors) == m]
  if (length(chems_mod)) {
    z_mat_sub <- z_mat_prod_all[chems_mod, , drop = FALSE]
    chem_id_lookup_mod <- data.frame(
      Label     = ifelse(nchar(rownames(z_mat_sub)) > 20,
                         paste0(substr(rownames(z_mat_sub), 1, 20), "…"),
                         rownames(z_mat_sub)),
      Chemical  = rownames(z_mat_sub),
      Module    = m,
      stringsAsFactors = FALSE
    )
    write.csv(chem_id_lookup_mod,
              file.path(out_dir_modsplit, sprintf("Heatmap_ChemID_lookup_%s.csv", m)),
              row.names = FALSE)
  }
}
