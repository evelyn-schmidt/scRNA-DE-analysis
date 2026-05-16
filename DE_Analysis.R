library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(presto)

merged <- readRDS('/cloud/project/data/single_cell_rna/backup_files/preprocessed_object.rds')

FeaturePlot(merged, features = 'Epcam')

DimPlot(merged, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  FeaturePlot(merged, features = 'Epcam') + 
  DimPlot(merged, group.by = 'immgen_singler_main')

VlnPlot(merged, group.by = 'seurat_clusters_res0.8', features = 'Epcam')

#set ident to seurat clusters metadata column and subset object to Epcam positive clusters
merged <- SetIdent(merged, value = 'seurat_clusters_res0.8')
merged_epithelial <- subset(merged, idents = c('9', '12'))

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(merged, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(merged_epithelial, group.by = 'seurat_clusters_res0.8', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(merged$seurat_clusters_res0.8)
table(merged_epithelial$seurat_clusters_res0.8)

#carry out DE analysis between both groups
merged_epithelial <- SetIdent(merged_epithelial, value = "seurat_clusters_res0.8")
epithelial_de <- FindMarkers(merged_epithelial, ident.1 = "9", ident.2 = "12", min.pct=0.25, logfc.threshold=0.1) #how cluster 9 changes wrt cluster 12


#restrict differentially expressed genes to those with an adjusted p-value less than 0.001 
epithelial_de_sig <- epithelial_de[epithelial_de$p_val_adj < 0.001,] 

#get the top 20 genes by fold change
epithelial_de_sig_top20 <- epithelial_de_sig %>%
  top_n(n = 20, wt = abs(avg_log2FC))

#get list of top 20 DE genes for ease
epithelial_de_sig_top20_genes <- rownames(epithelial_de_sig_top20)

#plot all 20 genes in violin plots
VlnPlot(merged_epithelial, features = epithelial_de_sig_top20_genes, 
        group.by = 'seurat_clusters_res0.8', ncol = 5, pt.size = 0)

#plot all 20 genes in UMAP plots
FeaturePlot(merged_epithelial, features = epithelial_de_sig_top20_genes, ncol = 5)

#plot all 20 genes in a DotPlot
DotPlot(merged_epithelial, features = epithelial_de_sig_top20_genes, 
        group.by = 'seurat_clusters_res0.8') + RotatedAxis()

#plot all differentially expressed genes in a volcano plot
EnhancedVolcano(epithelial_de,
                lab = rownames(epithelial_de),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Cluster9 wrt Cluster 12',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 0,
                colAlpha = 0.3)

#rerun FindMarkers
epithelial_de_gsea <- FindMarkers(merged_epithelial, ident.1 = "9", 
                                  ident.2 = "12", min.pct=0.25, logfc.threshold=0)
#save this table as a TSV file (first move index to first column)
epithelial_de_gsea <- tibble::rownames_to_column(epithelial_de_gsea, var = "gene")
write.table(x = epithelial_de_gsea, file = 'outdir_single_cell_rna/epithelial_de_gsea.tsv', sep='\t', row.names = FALSE)


# Aggregate counts for each sample and cluster combination
pb_epithelial <- AggregateExpression(merged_epithelial, assays = 'RNA', 
                                     return.seurat = T, group.by = c('orig.ident', 'seurat_clusters_res0.8'))
# See the first few rows of the log normalized data layer
print(head(pb_epithelial@assays$RNA$data))

#Change ident to seurat clusters
Idents(pb_epithelial) <- "seurat_clusters_res0.8"

# Use FindMarkers with DESeq2 as the test to compare cluster 9 wrt cluster 12
pb_epithelial_de <- FindMarkers(object = pb_epithelial, test.use = "DESeq2",
                                ident.1 = "9", ident.2 = "12")


# The DESeq2 analysis results in NAs in the pvalue columns for some cases 
pb_epithelial_de <- na.omit(pb_epithelial_de)

# Restrict differentially expressed genes to those with an adjusted p-value less than 0.001 
pb_epithelial_de_sig <- pb_epithelial_de[pb_epithelial_de$p_val_adj < 0.01,] 


# Compare significantly differentially expressed genes
pb_and_sc_genes <- intersect(rownames(epithelial_de_sig), rownames(pb_epithelial_de_sig))
only_sc_genes <- setdiff(rownames(epithelial_de_sig), rownames(pb_epithelial_de_sig))
only_pb_genes <- setdiff(rownames(pb_epithelial_de_sig), rownames(epithelial_de_sig))

print(paste0("Genes differentially expressed in both single-cell and pseudobulk: ", length(pb_and_sc_genes)))
print(paste0("Genes differentially expressed in single-cell but not pseudobulk: ", length(only_sc_genes)))
print(paste0("Genes differentially expressed in pseudobulk but not single-cell: ", length(only_pb_genes)))


#check all the annotated celltypes
unique(merged$immgen_singler_main)

#pick the ones that are related to T cells
t_celltypes_names <- c('T cells', 'NKT', 'Tgd')
merged <- SetIdent(merged, value = 'immgen_singler_main')
merged_tcells <- subset(merged, idents = t_celltypes_names)

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(merged, group.by = 'immgen_singler_main', label = TRUE) + 
  DimPlot(merged_tcells, group.by = 'immgen_singler_main', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(merged$immgen_singler_main)
table(merged_tcells$immgen_singler_main)

#we'll start by checking the possible names each replicate has.
unique(merged_tcells$orig.ident)

#there are 6 possible values, 3 replicates for the ICB treatment condition, and 3 for the ICBdT condition
#so we can combine "Rep1_ICB", "Rep3_ICB", "Rep5_ICB" to ICB, and "Rep1_ICBdT", "Rep3_ICBdT", "Rep5_ICBdT" to ICBdT. 
#first initialize a metadata column for experimental_condition
merged_tcells@meta.data$experimental_condition <- NA

#Now we can take all cells that are in each replicate-condition, 
#and assign them to the appropriate condition
merged_tcells@meta.data$experimental_condition[merged_tcells@meta.data$orig.ident %in% c("Rep1_ICB", "Rep3_ICB", "Rep5_ICB")] <- "ICB"
merged_tcells@meta.data$experimental_condition[merged_tcells@meta.data$orig.ident %in% c("Rep1_ICBdT", "Rep3_ICBdT", "Rep5_ICBdT")] <- "ICBdT"

#double check that the new column we generated makes sense 
#(each replicate should correspond to its experimental condition)
table(merged_tcells@meta.data$orig.ident, merged_tcells@meta.data$experimental_condition)

#carry out DE analysis between both groups
merged_tcells <- SetIdent(merged_tcells, value = "experimental_condition")
tcells_de <- FindMarkers(merged_tcells, ident.1 = "ICBdT", ident.2 = "ICB", min.pct=0.25)

#restrict differentially expressed genes to those with an adjusted p-value less than 0.001 
tcells_de_sig <- tcells_de[tcells_de$p_val_adj < 0.001,]

#find the top 5 most downregulated genes
tcells_de_sig %>%
  top_n(n = 5, wt = -avg_log2FC)
#find the top 5 most upregulated genes
tcells_de_sig %>%
  top_n(n = 5, wt = avg_log2FC)

# Subset object to CD8 T cells. Since we already showed how to subset cells using the clusters earlier, 
# This time we'll subset to CD8 T cells by selecting for cells with high expression of Cd8 genes and low expression of Cd4 genes using violin plots to find thresholds for filtering
VlnPlot(merged_tcells, features = c('Cd8a', 'Cd8b1', 'Cd4'))
merged_cd8tcells <- subset(merged_tcells, subset= Cd8b1 > 1 & Cd8a > 1 & Cd4 < 0.1)


#carry out DE analysis between both groups
merged_cd8tcells <- SetIdent(merged_cd8tcells, value = "experimental_condition")
cd8tcells_de <- FindMarkers(merged_cd8tcells, ident.1 = "ICBdT", ident.2 = "ICB", min.pct=0.25) #how ICBdT changes wrt ICB

#restrict differentially expressed genes to those with an adjusted p-value less than 0.001 
cd8tcells_de_sig <- cd8tcells_de[cd8tcells_de$p_val_adj < 0.001,]

#get the top 20 genes by fold change
cd8tcells_de_sig %>%
  top_n(n = 20, wt = abs(avg_log2FC)) -> cd8tcells_de_sig_top20

#get list of top 20 DE genes for ease
cd8tcells_de_sig_top20_genes <- rownames(cd8tcells_de_sig_top20)


#plot all 20 genes in violin plots
VlnPlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes, group.by = 'experimental_condition', ncol = 5, pt.size = 0)

DimPlot(merged_cd8tcells, group.by = 'experimental_condition')

#plot all 20 genes in UMAP plots
FeaturePlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes, ncol = 5)

#plot all 20 genes in a DotPlot
DotPlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes, group.by = 'experimental_condition') + RotatedAxis()

#plot all differentially expressed genes in a volcano plot
EnhancedVolcano(cd8tcells_de,
                lab = rownames(cd8tcells_de),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'ICBdT wrt ICB',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 5.0,
                colAlpha = 0.3)

#rerun FindMarkers
cd8tcells_de_gsea <- FindMarkers(merged_cd8tcells, ident.1 = "ICBdT", 
                                 ident.2 = "ICB", min.pct=0.25, logfc.threshold=0)
#save this table as a TSV file (first move index to first column)
cd8tcells_de_gsea <- tibble::rownames_to_column(cd8tcells_de_gsea, var = "gene")
write.table(x = cd8tcells_de_gsea, file = 'outdir_single_cell_rna/cd8tcells_de_gsea.tsv', sep='\t', row.names = FALSE)

# Aggregate counts for each sample and cluster combination
pb_cd8tcells <- AggregateExpression(merged_cd8tcells, assays = 'RNA', 
                                    return.seurat = T, group.by = c('orig.ident', 'experimental_condition'))

# See the first few rows of the log normalized data layer
print(head(pb_cd8tcells@assays$RNA$data))

# Change ident to seurat clusters
Idents(pb_cd8tcells) <- "experimental_condition"

# Use FindMarkers with DESeq2 as the test to compare cluster 9 wrt cluster 12
pb_cd8tcells_de <- FindMarkers(object = pb_cd8tcells, test.use = "DESeq2",
                               ident.1 = "ICBdT", ident.2 = "ICB")

# The DESeq2 analysis results in NAs in the pvalue columns for some cases 
pb_cd8tcells_de <- na.omit(pb_cd8tcells_de)

# Restrict differentially expressed genes to those with an adjusted p-value less than 0.001 
pb_cd8tcells_de_sig <- pb_cd8tcells_de[pb_cd8tcells_de$p_val_adj < 0.01,] 

# Compare significantly differentially expressed genes
pb_and_sc_genes <- intersect(rownames(cd8tcells_de_sig), rownames(pb_cd8tcells_de_sig))
only_sc_genes <- setdiff(rownames(cd8tcells_de_sig), rownames(pb_cd8tcells_de_sig))
only_pb_genes <- setdiff(rownames(pb_cd8tcells_de_sig), rownames(cd8tcells_de_sig))

print(paste0("Genes differentially expressed in both single-cell and pseudobulk: ", length(pb_and_sc_genes)))

print(paste0("Genes differentially expressed in single-cell but not pseudobulk: ", length(only_sc_genes)))
print(paste0("Genes differentially expressed in pseudobulk but not single-cell: ", length(only_pb_genes)))
