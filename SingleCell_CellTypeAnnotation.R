################################ Run Cell Annotation using a function ########################

singleCell_cellTypeAnnotation <- function(seurat_integrated) {
  # Cell type assignment
  # load gene set preparation function
  lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  # DB file
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue = "Immune system" # e.g. Immune   system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  
  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seurat_integrated[["integrated"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # Resolution 1.4(seurat_integrated@meta.data$integrated_snn_res.1.4)
  cL_results = do.call("rbind", lapply(unique(seurat_integrated@meta.data$integrated_snn_res.1.4), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$integrated_snn_res.1.4==cl, ])]),   decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_integrated@meta.data$integrated_snn_res.1.4==cl)), 5)
  }))
  
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  write.csv(sctype_scores, "sctype_scores_final_1.4.csv")
  write.csv(cL_results, "res1.4_bubble_cluster_scores.csv")
  
  # Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
  n_cells <- FetchData(seurat_integrated, 
                       vars = c("ident", "sample")) %>%
    group_by(sample) %>%
    dplyr::count(ident) %>% 
    spread(ident, n) 
  
  write.csv(n_cells, "res1.4_cells_per_cluster.csv")
}

# Setup the current working directory
setwd("L:/BIOINFORMATICS_SHARED/Yan/Harry_Kim/Y276_Res1.4_Seurat_output.tar/Y276_Res1.4_Seurat_output/Seurat_output/data/results")

# Read Seurat object RDS file.
seurat_integrated_Y276 <- readRDS("Y276_final_integrated_seurat.rds")

# Call the singleCell_cellTypeAnnotation function
singleCell_cellTypeAnnotation(seurat_integrated_Y276)

