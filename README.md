# Transcriptomic and Functional Responses of Human Airway Cells to Vaped ∆8-THC and its Oxidation Product ∆8-THCQ

This code was generated to support the manuscript entitled, 'Transcriptomic and Functional Responses of Human Airway Cells to Vaped ∆8-THC and its Oxidation Product ∆8-THCQ' currently accepted for publication at Scientific Reports. doi: https://doi.org/10.1038/s41598-026-53356-z. 

> This project investigates the respiratory health risks of commercial Δ8-THC vape products. We identify Δ8-THC quinone (Δ8-THCQ, HU-336) as a major constituent of Δ8-THC distillates and disposables, with concentrations rising ~3.67-fold to millimolar levels after vaping. Using bronchial epithelial cell models and the UNC Vaping Product Exposure System (VaPES), we show that Δ8-THCQ and Δ8-THC aerosols activate stress, inflammatory, and fibrosis-linked pathways while impairing motile cilia function. Together, these findings suggest repeated Δ8-THC vaping may disrupt mucociliary defense and raise concern for chronic airway injury. 


# D8_D8Q_16HBE_RNA-SeqAnalysis

- Using DESeq2 to determine differentially expressed genes after varying Δ8-THC products (control, Δ8-THC, and Δ8-THCQ) in 16HBE bronchial epithelial cells at 12 and 24 hour exposures
- Creating PCA plots to visualize sample clustering across treatment groups
- Creating heatmaps of genes that were significantly different across conditions
- Creating Venn/Euler diagrams to determine genes significantly altered across multiple exposure comparisons
- Generating volcano plots of differentially expressed genes with top gene labels


# VapedD8_HBEC_RNA-SeqAnalysis

- Using limma/voom to determine differentially expressed genes in primary human bronchial epithelial cells exposed to aerosols from six commercial Δ8-THC products (Palmetto Distillate, NYSW Distillate, Palmetto Disposable, Stiizy Disposable, ST Strawberry Juice, ST Green Apple Juice) vs. PG/VG control
- Accounting for technical replicates with duplicateCorrelation and correcting for donor batch effects with removeBatchEffect
- Creating MDS plots before and after donor batch correction with k-means cluster ellipses
- Creating volcano plots, heatmaps, and Venn/Euler diagrams comparing significant genes across product types
- Running linear dose-response modeling against measured Δ8-THC and Δ8-THCQ concentrations to identify dose-responsive genes, with top-gene dose-response curves
- Comparing dose-responsive genes from vaped Δ8-THC exposures to Δ8-THCQ-responsive genes from apical 16HBE exposures via Euler diagrams


# Vaped_D8_WCCNA

- Using Weighted Chemical Correlation Network Analysis (WCCNA) to correlate the chemical composition of Δ8-THC vaped condensates (GC-MS/LC-MS data) with RNA-seq data from primary differentiated human bronchial epithelial cells exposed to Δ8-THC aerosols
- Filtering and normalizing RNA-seq data with DESeq2 (collapsing technical replicates, retaining genes significantly altered in ≥3 exposures vs. control)
- Constructing a chemistry co-expression network, detecting chemical modules, and computing module eigengenes
- Correlating module eigengenes with RNA-seq traits to identify chemical clusters driving specific airway transcriptional responses
- Generating module-trait correlation heatmaps, chemical abundance heatmaps grouped by module and product, and per-module composite heatmaps pairing compound abundances with top gene correlations
