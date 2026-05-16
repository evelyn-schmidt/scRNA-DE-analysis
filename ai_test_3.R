# ============================================================
# CD8 T Cell DE Analysis
# Method: Wilcoxon rank-sum
# Condition column: experimental_condition
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# ── 0. Load your Seurat object ───────────────────────────────
# Replace this line with however you load your object, e.g.:
seurat_obj <- readRDS("/cloud/project/data/single_cell_rna/backup_files/preprocessed_object.rds")
seurat_obj@meta.data$experimental_condition[seurat_obj@meta.data$orig.ident %in% c("Rep1_ICB", "Rep3_ICB", "Rep5_ICB")] <- "ICB"
seurat_obj@meta.data$experimental_condition[seurat_obj@meta.data$orig.ident %in% c("Rep1_ICBdT", "Rep3_ICBdT", "Rep5_ICBdT")] <- "ICBdT"


### TRY IT WITH SETTING THE THRESHOLD
# 1. Visual check to see where your CD8 cells are hiding
FeaturePlot(seurat_obj, features = c("Cd3e", "Cd8a", "Cd8b"))


# 1. First, pull out ALL T cells using a broad marker like CD3E
# Adjust the expression threshold (e.g., > 0.5) based on your data's normalization
all_t_cells <- subset(seurat_obj, subset = Cd3e > 0.5)

# 2. Plot CD8A vs CD4 to see the separation
FeatureScatter(all_t_cells, feature1 = "Cd8a", feature2 = "Cd4")

# 3. Strictly subset the cells that are CD8+ and CD4-
cd8_cells <- subset(all_t_cells, subset = Cd8a > 0.5  & Cd4 == 0)

# Ensure you are analyzing the CD8 subset
DefaultAssay(cd8_cells) <- "RNA"

# Set the active identity to your treatment column
# (Replace "treatment_group" with your actual metadata column name, e.g., "Condition")
Idents(cd8_cells) <- "experimental_condition"

# Run DE using MAST (Make sure library(MAST) is installed via BiocManager)
de_results <- FindMarkers(
  object = cd8_cells,
  ident.1 = "ICBdT",      # Name of your treatment group
  ident.2 = "ICB",      # Name of your control group
  test.use = "wilcox",
)

library(EnhancedVolcano)

EnhancedVolcano(
  de_results,
  lab = rownames(de_results),
  x = 'avg_log2FC',
  y = 'p_val_adj', # Using adjusted p-value for statistical stringency
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 2.0,
  labSize = 4.0,
  title = 'CD8+ T Cells: Treated vs Control',
  subtitle = 'Cell-level DE analysis',
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 3.0,
  widthConnectors = 0.5
)

# 1. Filter for statistically significant genes first (e.g., adjusted p-value < 0.05)
strict_de <- de_results[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1, ]

# 2. Get the Top 5 UP-regulated genes (highest positive log2FC)
head(strict_de[order(-strict_de$avg_log2FC), ], 5)
head(rownames(strict_de[order(-strict_de$avg_log2FC), ]), 5)

# 3. Get the Top 5 DOWN-regulated genes (lowest negative log2FC)
head(strict_de[order(strict_de$avg_log2FC), ], 5)
head(rownames(strict_de[order(strict_de$avg_log2FC), ]), 5)



#### Try with cell type labeling
# Install SingleR and reference datasets if needed
BiocManager::install(c("SingleR", "celldex"))
library(SingleR)
library(celldex)

# 1. Load a high-quality immune reference dataset
ref <- MonacoImmuneData() # Excellent for detailed T cell subsets

# 2. Convert Seurat object to SingleCellExperiment format for SingleR
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seurat_obj)

# 3. Run the annotation
predictions <- SingleR(test = sce, ref = ref, labels = ref$label.fine)

# 4. Add the predictions back to your Seurat metadata
your_seurat_object$singleR_labels <- predictions$pruned.labels

# 5. Subset specifically for the CD8 categories identified by SingleR
cd8_cells <- subset(your_seurat_object, subset = singleR_labels %in% c("CD8+ T cells", "Naive CD8+ T cells", "Effector memory CD8+ T cells"))

de_results <- FindMarkers(
  object = cd8_cells,
  ident.1 = "ICBdT",      # Name of your treatment group
  ident.2 = "ICB",      # Name of your control group
  test.use = "wilcox",
)

EnhancedVolcano(
  de_results,
  lab = rownames(de_results),
  x = 'avg_log2FC',
  y = 'p_val_adj', # Using adjusted p-value for statistical stringency
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 2.0,
  labSize = 4.0,
  title = 'CD8+ T Cells: Treated vs Control',
  subtitle = 'Cell-level DE analysis',
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 3.0,
  widthConnectors = 0.5
)

# 1. Filter for statistically significant genes first (e.g., adjusted p-value < 0.05)
strict_de <- de_results[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1, ]

# 2. Get the Top 5 UP-regulated genes (highest positive log2FC)
head(strict_de[order(-strict_de$avg_log2FC), ], 5)
head(rownames(strict_de[order(-strict_de$avg_log2FC), ]), 5)

# 3. Get the Top 5 DOWN-regulated genes (lowest negative log2FC)
head(strict_de[order(strict_de$avg_log2FC), ], 5)
head(rownames(strict_de[order(strict_de$avg_log2FC), ]), 5)


