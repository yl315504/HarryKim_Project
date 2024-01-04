# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_05_dge.html#Pseudobulk
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RCurl))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(edgeR))

# Create the directory "DEG" for the first time.
dir.create("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG")

# The function to do Pseudobulk differential expression analysis
DEG <- function(seurat_object, cell_type){
  sub_data <- seurat_object
  DGE_DATA <- sub_data@assays$RNA@counts
  
  # Create a sparse model matrix for pseudobulk calculation
  mm <- Matrix::sparse.model.matrix(~0 + sub_data$orig.ident)
  
  # Calculate the pseudobulk expression by multiplying counts and model matrix
  pseudobulk <- DGE_DATA %*% mm
  
  # Define labels for bulk samples
  # Modify the labels based on your study design
  bulk.labels = c("left", "right", "left", "right", "left", "right", "left", "right")
  
  # Create a DGEList object with pseudobulk expression and bulk labels
  dge.list <- DGEList(counts = pseudobulk, group = factor(bulk.labels))
  
  # Filter out lowly expressed genes
  keep <- filterByExpr(dge.list)
  dge.list <- dge.list[keep, , keep.lib.sizes = FALSE]
  
  # Normalize expression values
  dge.list <- calcNormFactors(dge.list)
  
  # Define the design matrix for the GLM model
  design = model.matrix(~bulk.labels)
  
  # Estimate dispersion
  dge.list <- estimateDisp(dge.list, design)
  
  # Fit a negative binomial GLM model
  fit <- glmQLFit(dge.list, design)
  
  # Perform a likelihood ratio test
  qlf <- glmQLFTest(fit, coef = 2)
  
  # Extract all and top 500 differential expressed genes
  res.edgeR <- topTags(qlf, n = Inf)$table
  # res.edgeR.top500 <- topTags(qlf, n = 500)$table
  
  # Extract and process the results from the likelihood ratio test
  res.edgeR$gene = rownames(res.edgeR)
  res.edgeR$dir = ifelse(res.edgeR$logFC > 0, "left", "right")
  res.edgeR.top500$gene = rownames(res.edgeR.top500)
  
  # Select the top differential expressed genes based on direction and p-value
  res.edgeR %>%
    group_by(dir) %>%
    top_n(-10, PValue) %>%
    arrange(dir) -> top.edgeR
  
  # Create a DotPlot using Seurat's plotting functions
  out_dir <- "D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\"
  
  png(paste0(out_dir, cell_type,"_DEG.png"), width = 12, height = 12, units = "in", res=300)
  output_plot <- DotPlot(
    seurat_object,
    features = as.character(unique(top.edgeR$gene)),
    group.by = "orig.ident",
    assay = "RNA"
  ) + coord_flip() + ggtitle(paste0(cell_type," pseudobulk")) + RotatedAxis()
  print(output_plot)
  dev.off()
  
  tmp1 <- paste0(out_dir, cell_type, "_DEG.txt")
  write.table(res.edgeR %>% select(logFC, PValue, FDR, gene), tmp1, sep = "\t", row.names = FALSE)
  
  tmp2 <- paste0(out_dir, cell_type, "_DEG_top500.txt")
  write.table(res.edgeR.top500 %>% select(logFC, PValue, FDR, gene), tmp2, sep = "\t", row.names = FALSE)
}

setwd("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results")
seurat_integrated_CellAnnotation_FourSamples <- readRDS("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\seurat_integrated_CellAnnotation_FourSamples.rds")

##################################### How many cells are in each cluster ###############################################
# Idents(seurat_integrated_CellAnnotation_FourSamples) <- "cell_type"
CellNum_per_CellType <- table(Idents(seurat_integrated_CellAnnotation_FourSamples), seurat_integrated_CellAnnotation_FourSamples$orig.ident)
dir_name <- "D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\cell_number_per_type.txt"
write.table(CellNum_per_CellType, dir_name, sep = "\t", row.names = TRUE)

##################################### DGE between left and right in all clusters (cell types) #########################
# Two ways to get subsets
# seurat_Macrophages <- subset(seurat_integrated_CellAnnotation_FourSamples, subset = cell_type == "Macrophages")
# DEG(seurat_Macrophages, "Macrophages")
# 
# seurat_T_cells <- subset(seurat_integrated_CellAnnotation_FourSamples, idents = c("Naive CD4+ T cells", "Effector CD8+ T cells","Memory CD4+ T cells", "CD8+ NKT-like cells" ))
# DEG(seurat_T_cells, "T Cells")

# cell_list has all cell types. cell_list is used in the for loop later
cell_list <- unique(seurat_integrated_CellAnnotation_FourSamples@active.ident)

# Run the DEG function on all cell types using a for loop
for (cell in cell_list){
  seurat_subset <- subset(seurat_integrated_CellAnnotation_FourSamples, idents = c(cell))
  DEG(seurat_subset, cell)
}


# Gene Set Enrichment Analysis (GSEA): Compute the score for each gene for Macrophages
# Macrophages <- read.delim("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\Macrophages_DEG.txt", header = TRUE, sep = "\t")
# Macrophages$rankMetric <- sign(Macrophages$logFC)*(-log10(Macrophages$PValue))
# write.table(Macrophages[, c("gene", "rankMetric")], file = "D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\GSEA_score_Macrophages.rnk", sep = "\t",
#             quote = FALSE, row.names = FALSE, col.names = FALSE) 
 
# Gene Set Enrichment Analysis (GSEA): Compute the score for each gene for all cell types using a for loop
for (cell in cell_list){
  CellType <- cell
  cell <- read.delim(paste0("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\", CellType, "_DEG.txt"), header = TRUE, sep = "\t")
  cell$rankMetric <- sign(cell$logFC)*(-log10(cell$PValue))
  write.table(cell[, c("gene", "rankMetric")], file = paste0("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\pathway_results\\GSEA_score_", CellType, ".rnk"), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE) 
}

# Gene Set Enrichment Analysis (GSEA): Compute the score for each gene for all T cells
T_cells <- read.delim("D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\T_cells_DEG.txt", header = TRUE, sep = "\t")
T_cells$rankMetric <- sign(T_cells$logFC)*(-log10(T_cells$PValue))
write.table(T_cells[, c("gene", "rankMetric")], file = "D:\\Yan\\Harry_Kim\\four_samples\\Seurat_output\\data\\results\\DEG\\pathway_results\\GSEA_score_T_cells.rnk", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



