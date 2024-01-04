# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RCurl))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(AnnotationHub))

in_dir = "L:/BIOINFORMATICS_SHARED/Yan/Harry_Kim/all_samples"
SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
out_prefix = "all_samples"
org = "human"

setwd(in_dir)
out_dir <- paste0(in_dir,"/Seurat_output/data/results/")
# dir.create(out_dir, recursive = TRUE)
# dir.create(paste0(out_dir,"QC/"))
# 
# SamplesSheet <- read.csv(SamplesSheet_file,header=TRUE)
# 
# files <- list.files(path = paste0(in_dir,"/Seurat_input_data/"), pattern = "*.h5", full.names = T)

# paste0("Input Directory: ",in_dir)
# paste0("Out prefix: ",out_prefix)
# paste0("SamplesSheet file: ",SamplesSheet_file)
# paste0("Output Directory: ",out_dir)
# paste0("Output QC Directory: ",out_dir,"QC/")
# 
# paste0("SamplesSheet input: \n")
# SamplesSheet
# paste0("Input Hdf5r files: \n")
# files

# Read filtered after QC seurat object
filtered_seurat <- readRDS(file=paste0(out_dir,"seurat_filtered.rds"),ref=NULL)

# Single-cell RNA-seq analysis - clustering analysis

# if (org == "mouse") {
#   m.cc.genes <- readRDS("/home1/08270/yanliu74/Seurat_scRNASeq/mouse_cell_cycle_genes.rds")
#   s_genes <- m.cc.genes$s.genes
#   g2m_genes <- m.cc.genes$g2m.genes
# } else if (org == "human") {
#   # cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
#   # cell_cycle_genes <- read.csv(text = Homo_sapiens.csv)
#   cell_cycle_genes <- read.csv("L:/BIOINFORMATICS_SHARED/Yan/Harry_Kim/all_samples/Homo_sapiens.csv")
#   # Connect to AnnotationHub
#   ah <- AnnotationHub()
#   
#   # Access the Ensembl database for organism
#   ahDb <- query(ah, 
#                 pattern = c("Homo sapiens", "EnsDb"), 
#                 ignore.case = TRUE)
#   
#   # Acquire the latest annotation files
#   id <- ahDb %>%
#     mcols() %>%
#     rownames() %>%
#     tail(n = 1)
#   
#   # Download the appropriate Ensembldb database
#   edb <- ah[[id]]
#   
#   # Extract gene-level information from database
#   annotations <- genes(edb, 
#                        return.type = "data.frame")
#   
#   # Select annotations of interest
#   annotations <- annotations %>%
#     dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
#   
#   # Get gene names for Ensembl IDs for each gene
#   cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
#   
#   # Acquire the S phase genes
#   s_genes <- cell_cycle_markers %>%
#     dplyr::filter(phase == "S") %>%
#     pull("gene_name")
#   
#   # Acquire the G2M phase genes        
#   g2m_genes <- cell_cycle_markers %>%
#     dplyr::filter(phase == "G2/M") %>%
#     pull("gene_name")
# } else {
#   print("Organism other than human or mouse")
# }

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

# Score cells for cell cycle
# seurat_phase <- CellCycleScoring(seurat_phase, 
#                                  g2m.features = g2m_genes, 
#                                  s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
head(seurat_phase@meta.data)     

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
# png(paste0(out_dir,"QC/After_filter/",out_prefix,"_cellcyclephase.png"), width = 12, height = 12, units = 'in', res = 300)    
# DimPlot(seurat_phase,
#         reduction = "pca",
#         group.by= "Phase",
#         split.by = "Phase")
dev.off()

# options(future.globals.maxSize = 4000 * 1024^2)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

#split_seurat <- split_seurat[c("ctrl", "stim")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Check which assays are stored in objects
split_seurat[[1]]@assays

# Save individual seurat object
lapply(1:length(split_seurat), function(i) saveRDS(split_seurat[[i]], paste0(out_dir,"/",names(split_seurat)[i],"_seurat.rds")))

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
k.filter <- min(200, min(sapply(split_seurat, ncol)))

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", k.filter = k.filter,
                                        anchor.features = integ_features)

memory.limit(56000)									
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

