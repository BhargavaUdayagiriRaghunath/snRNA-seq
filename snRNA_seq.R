
#install.packages("Seurat")
library(Seurat)

# Loading count matrix data
data_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138852&format=file&file=GSE138852%5Fcounts%2Ecsv%2Egz"
download.file(data_url, "GSE138852_counts.csv.gz")
data <- read.table("GSE138852_counts.csv", sep = ",", header = TRUE, row.names = 1)

# Creating a Seurat object
seurat_obj <- CreateSeuratObject(counts = data)

#VlnPlot(seurat_obj, features = c("DPP10", "CEMIP", "SLC1A2"), ncol = 3)
all.genes <- rownames(seurat_obj)

# Preprocessing and Quality Control
  # SCTransform used to normalize the data 
  # Principal Component Analysis is run to lessen the noise and reduce dimensionality
  # Constructing Shared nearest neighbors graph
  # Forming clusters from SNN graph

seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.4)  # Resolution set to 1.4, since nuclei count > 5000

# Viewing columns associated with the metadata in the seurat object
View(seurat_obj@meta.data)

# Visualization by utilization of UMAP and DimPlot to visualize
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE)

# Differential Expression Analysis

clusters <- unique(seurat_obj$seurat_clusters)  # factoring all different clusters

degs_by_cluster <- list() # Creating empty list to store list of each cluster

  # Loop iterating over all the 25 clusters
for (cluster_id in clusters) {
  deg_markers <- FindMarkers(seurat_obj, 
                          test.use = "MAST", 
                          ident.1 = cluster_id,
                          logfc.threshold = 0.25, 
                          min.pct = 0.1, assay = "RNA")
  degs_by_cluster[[as.character(cluster_id)]] <- deg_markers
}
head(seurat_obj)

# Print each of the differentially expressed genes in each cluster
for (cluster_id in clusters) {
  cluster_degs <- degs_by_cluster[[as.character(cluster_id)]]
  cat("Cluster", cluster_id, "DEGs:\n")
  print(cluster_degs)
}



#saveRDS(seurat_obj, "seurat_object.rds")
