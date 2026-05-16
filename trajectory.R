library(Seurat)
library(CytoTRACE)
library(monocle3)
library(ggplot2)

devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") #installing
library(CytoTRACE2) #loading

rep135 <- readRDS("data/single_cell_rna/backup_files/preprocessed_object.rds")

DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE) + 
  DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) 

Epithelial_cells = rep135$immgen_singler_main =="Epithelial cells"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Epithelial cells")
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)

FeaturePlot(object = rep135, features = c("Cd44", "Krt14", "Krt5", "Krt16", "Krt6a"))

cell_type_Basal_marker_gene_list <- list(c("Cd44", "Krt14", "Krt5", "Krt16", "Krt6a"))
rep135 <- AddModuleScore(object = rep135, features = cell_type_Basal_marker_gene_list, name = "cell_type_Basal_score") 
FeaturePlot(object = rep135, features = "cell_type_Basal_score1")

cell_type_Luminal_marker_gene_list <- list(c("Cd24a", "Erbb2", "Erbb3", "Foxa1", "Gata3", "Gpx2", "Krt18", "Krt19", "Krt7", "Krt8", "Upk1a"))
rep135 <- AddModuleScore(object = rep135, features = cell_type_Luminal_marker_gene_list, name = "cell_type_Luminal_score")
FeaturePlot(object = rep135, features = "cell_type_Luminal_score1")

### Subsetting dataset epithelial
rep135 <- SetIdent(rep135, value = 'seurat_clusters_res0.8')
rep135_epithelial <- subset(rep135, idents = c('9', '12')) # 1750

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_epithelial, group.by = 'seurat_clusters_res0.8', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(rep135$seurat_clusters_res0.8)
table(rep135_epithelial$seurat_clusters_res0.8)

rep135_epithelial_expression <- data.frame(GetAssayData(object = rep135_epithelial, layer = "data"))

rep135_epithelial_cytotrace_scores <- CytoTRACE(rep135_epithelial_expression, ncores = 1)

rep135_epithelial_cytotrace_transposed <- as.data.frame(rep135_epithelial_cytotrace_scores$CytoTRACE) 
names(rep135_epithelial_cytotrace_transposed) <- "cytotrace_scores"

head(rep135_epithelial_cytotrace_transposed)

# fix the barcode formatting 
rownames(rep135_epithelial_cytotrace_transposed) <- sub("\\.", "-", rownames(rep135_epithelial_cytotrace_transposed))
rownames(rep135_epithelial_cytotrace_transposed)

rep135_epithelial <- AddMetaData(rep135_epithelial, rep135_epithelial_cytotrace_transposed %>% select("cytotrace_scores"))

rep135_epithelial[['differentiation_scores']] <- 1 - rep135_epithelial[['cytotrace_scores']] # Let's also reverse out CytoTRACE scores so that high means more differentiated and low means less differentiated

FeaturePlot(object = rep135_epithelial, features = c("cell_type_Basal_score1", "cell_type_Luminal_score1", "cytotrace_scores", "differentiation_scores"))

#### CytoTrace2

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
#rep135_epithelial_expression <- data.frame(GetAssayData(object = rep135_epithelial, layer = "counts"))
cytotrace2_result <- cytotrace2(rep135_epithelial, is_seurat = TRUE, slot_type = "counts", species = 'mouse')
cytotrace2_result

# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = rep135_epithelial@meta.data$seurat_clusters_res0.8) %>% set_rownames(., colnames(rep135_epithelial))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)

plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Boxplot_byPheno
plots$Phenotype_UMAP

### Just with the liminal Cells

DimPlot(rep135_epithelial) # cluster 10 is our luminal cells

rep135_luminal <- subset(rep135_epithelial, idents = c('9')) # 863 cells

rep135_luminal_expression <- data.frame(GetAssayData(object = rep135_luminal, layer = "data"))

rep135_luminal_cytotrace_scores <- CytoTRACE(rep135_luminal_expression, ncores = 1)

rep135_luminal_cytotrace_transposed <- as.data.frame(rep135_luminal_cytotrace_scores$CytoTRACE) 
names(rep135_luminal_cytotrace_transposed) <- "cytotrace_scores"

head(rep135_luminal_cytotrace_transposed)

rownames(rep135_luminal_cytotrace_transposed) <- sub("\\.", "-", rownames(rep135_luminal_cytotrace_transposed))
rownames(rep135_luminal_cytotrace_transposed)

rep135_luminal <- AddMetaData(rep135_luminal, rep135_luminal_cytotrace_transposed)

rep135_luminal[['cytotrace_scores_luminal']] <- 1 - rep135_luminal[['cytotrace_scores']]
rep135_luminal[['differentiation_scores_luminal']] <- 1 - rep135_luminal[['differentiation_scores']]

# compare all epithelial cells CytoTRACE scores to the luminal-only CytoTRACE
(FeaturePlot(object = rep135_luminal, features = c("cytotrace_scores_luminal")) +
    ggtitle("Luminal Cells CytoTRACE Scores")) +
  (FeaturePlot(object = rep135_luminal, features = c("differentiation_scores_luminal")) +
     ggtitle("Luminal Cells Differentiation Scores"))

# View the S-phase genes, G2/M-phase genes, and the Phase to see if that explains the differentiation score
FeaturePlot(object = rep135_luminal, features = c("S.Score", "G2M.Score")) + 
  DimPlot(rep135_luminal, group.by = "Phase")

FeaturePlot(object = rep135_luminal, features = c("S.Score", "G2M.Score", "differentiation_scores"))

rep135_cds <- SeuratWrappers::as.cell_data_set(rep135)
rep135_cds <- cluster_cells(rep135_cds)
plot_cells(rep135_cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

rep135_cds <- learn_graph(rep135_cds, use_partition = FALSE) # graph learned across all partitions

rep135_cds <- order_cells(rep135_cds)

plot_cells(rep135_cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE, cell_size = 1)

rep135_cytotrace <- read.table("rep135_cytotrace_scores.tsv", sep="\t") 
head(rep135_cytotrace)

rep135_cytotrace_transposed <- t(rep135_cytotrace)
colnames(rep135_cytotrace_transposed) <- "CytoTRACE"
head(rep135_cytotrace_transposed)

rownames(rep135_cytotrace_transposed) <- sub("\\.", "-", rownames(rep135_cytotrace_transposed))
rownames(rep135_cytotrace_transposed)

rep135 <- AddMetaData(rep135, rep135_cytotrace_transposed)
rep135[["differentiation_score"]] <- 1 - rep135[["CytoTRACE"]]

plot_cells(rep135_cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE, cell_size = 1) + 
  (FeaturePlot(rep135, features = 'differentiation_score') + scale_color_viridis(option = 'magma', discrete = FALSE)) +
  DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE)


highlight = rep135$immgen_singler_main =="Macrophages"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Macrophages")
# Plot the UMAP
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)


highlight = rep135$immgen_singler_main =="Monocytes"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Monocytes")
# Plot the UMAP
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)


# grab all cells that are macrophages and monocytes
Idents(rep135) <- "immgen_singler_main" 
rep135_macro_mono_cells <- subset(rep135, idents = c("Macrophages", "Monocytes"), invert = FALSE) # 1092

DimPlot(rep135_macro_mono_cells, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) 

# grab all cells that are macrophages and monocytes, we can subset by clusters 6 and 14 which seem to contain 
Idents(rep135) <- "seurat_clusters_res0.8" 
rep135_macro_mono_cells <- subset(rep135, idents = c(6, 14), invert = FALSE) # 1350

DimPlot(rep135_macro_mono_cells, group.by = 'seurat_clusters', label = TRUE) + 
  DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) 

DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_fine', label = TRUE)


# create a data frame with the counts from our subsetted obect
rep135_macro_mono_cells_expression <- data.frame(GetAssayData(rep135_macro_mono_cells, layer = "data")) 

# pass that dataframe to the CytoTRACE function
rep135_macro_mono_cells_cytotrace_scores <- CytoTRACE(rep135_macro_mono_cells_expression, ncores = 4)

# Create a dataframe out of the CytoTRACE scores
rep135_macro_mono_cells_cytotrace_scores_df <- as.data.frame(rep135_macro_mono_cells_cytotrace_scores$CytoTRACE)

# Make the rownames of the cytotraace scores function the cell barcodes and rename the CytoTRACE scores column approproately
rownames(rep135_macro_mono_cells_cytotrace_scores_df) <- sub("\\.", "-", rownames(rep135_macro_mono_cells_cytotrace_scores_df))