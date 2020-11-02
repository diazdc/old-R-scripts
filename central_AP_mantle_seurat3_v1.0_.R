# 1.0 Subsetted cell types from Nmast_sub_anchored_seurat3_v1.0
{
library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)

# Clear workspace
rm(list = ls())

script_name <- paste0("central_AP_mantle_seurat3_v1.0")

dir.create(paste0("/n/projects/ddiaz/Analysis",
  "/Scripts/sb2191-regen/Seurat3-final/", script_name,"_figures/"))

figurePath <- function(filename){paste0("/n/projects/ddiaz/Analysis",
  "/Scripts/sb2191-regen/Seurat3-final/", script_name, "_figures/", filename)}

dataPath <- function(filename){paste0("/n/projects/",
  "ddiaz/Analysis/Data/sb2191-regen/", filename)}

gene_info <- read.delim(paste0("/n/projects/ddiaz/Analysis/",
  "Data/gene-lists/Danio_Features_unique_Ens91_v2.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

trt_colors <- c("green3", "gold", "darkorange", "red", "magenta",
  "mediumpurple1", "lightseagreen", "deepskyblue", "blue")
}


# ============================================================ Reading in Data


seurat_obj <- readRDS(paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "SeurObj_anchored_cell_type_update_nmast_only_additional_timepoints",
  "_seurat3_v1.0_.RDS"))

seurat_obj@meta.data$cell.type.ident <- plyr::revalue(
  seurat_obj@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))

DefaultAssay(seurat_obj) <- "RNA"

meta <- seurat_obj@meta.data
types <- levels(meta$cell.type.ident)
cell_names <- types
obj_list <- list()[1:8]

for (i in seq_along(types)){
  obj_list[[i]] <- seurat_obj[,meta$cell.type.ident == types[i]]
}

obj_list <- setNames(obj_list, types)


# ==================================================================== Gene QC


for (i in seq_along(types)) {
  obj_list[[i]] <- NormalizeData(obj_list[[i]],
  normalization.method = "LogNormalize", scale.factor = 10000)

  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]],
  selection.method = "vst", nfeatures = 2000)

  top10 <- head(VariableFeatures(obj_list[[i]]), 10)
  plot1 <- VariableFeaturePlot(obj_list[[i]])
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

  pdf(figurePath(paste0(cell_names[1], "_gene_QC.pdf")),
    width = 20, height = 10)
  CombinePlots(plots = list(plot1, plot2))
  dev.off()

  all_genes <- rownames(obj_list[[i]])
  obj_list[[i]] <- ScaleData(obj_list[[i]], features = all_genes)
  print(obj_list[[i]])
}


# ============================================================== Dim reduction


for (i in seq_along(types)) {
  obj_list[[i]] <- RunPCA(obj_list[[i]],
    features = VariableFeatures(object = obj_list[[i]]))

  pc_dims <- 1:50
  if (FALSE) { 
      pdf(figurePath(paste0(cell_names[i], "_PC_elbowPlot.pdf")))
      ElbowPlot(obj_list[[i]], ndims = 50)
      dev.off()
    
      png(figurePath(paste0(cell_names[i], "_PC_hmap1-15.png")),
        width = 20, height = 50, units = "in", res = 200)
      DimHeatmap(obj_list[[i]], dims = 1:15, nfeatures = 100,cells = 1000)
      dev.off()
    
      png(figurePath(paste0(cell_names[i], "_PC_hmap16-30.png")),
        width = 20, height = 50, units = "in", res = 200)
      DimHeatmap(obj_list[[i]], dims = 16:30, nfeatures = 100, cells = 1000)
      dev.off()
      print(obj_list[[i]])
  }
}


# ======================================================== Clustering and UMAP


dims <- 1:8
res = 0.6

for (i in seq_along(types)) {
  obj_list[[i]] <- FindNeighbors(obj_list[[i]], dims = dims, k.param = 20)
  obj_list[[i]] <- FindClusters(obj_list[[i]], resolution = res)
  obj_list[[i]] <- BuildClusterTree(obj_list[[i]],
    reorder = TRUE, reorder.numeric = TRUE)

  obj_list[[i]] <- RunUMAP(obj_list[[i]],
    dims = dims, min.dist = 0.3, spread = 1)
}


# ==== Plot by cluster
for (i in seq_along(types)) {
  umap_clusters <- DimPlot(obj_list[[i]], reduction = "umap", pt.size = 1.80,
    label = FALSE, label.size = 4)

umap_clusters <- cleanUMAP(umap_clusters)

  try({png(figurePath(paste0(types[i], "/",
      cell_names[i], "_UMAP_clusters_", max(dims), "_res_", res,
      "_5.png")), width = 10, height = 8, units = "in", res = 300)
    print(umap_clusters)
    dev.off()})
}


# ==== Plot by treatment
for (i in seq_along(types)) {
  umap_dataset <- DimPlot(obj_list[[i]], reduction = "umap", pt.size = 1.80,
    label = FALSE, label.size = 4, group.by = "data.set", cols = trt_colors)

  umap_dataset <- cleanUMAP(umap_dataset)

  try({png(figurePath(paste0(types[i],"/",
    cell_names[i], "_UMAP_dataset_", max(dims),"_5.png")),
    width = 10, height = 8, units = "in", res = 300)
    print(umap_dataset)
    dev.off()})
}



# ============================================================= Save/load data


if(TRUE) {
  saveRDS(obj_list,
    dataPath(paste0("CombinedObj_", script_name,"_.RDS")))
  obj_list <- readRDS(
    dataPath(paste0("CombinedObj_", script_name,"_.RDS")))
  types <- levels(obj_list[[1]]@meta.data$cell.type.ident)
  cell_names <- types
}

colnames(obj_list[[1]]@meta.data)

for (i in seq_along(types)){
  str(obj_list[[i]]@assays$RNA@counts)
}

table(seurat_obj@meta.data$data.set)
stop("break section 2")


# ========================================================= Additional figures
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

types <- levels(obj_list[[1]]@meta.data$cell.type.ident)
for (i in 1:4) {
  dir.create(figurePath(types[i]), showWarnings = FALSE)
}

for (i in seq_along(obj_list)) {
  DefaultAssay(obj_list[[i]]) <- "RNA"
}

# ==== mclapply for Find Markers
marker_table_list <- list()
for (i in seq_along(obj_list)) {
  marker_table <- mcFindMarkers(obj_list[[i]],
    save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "")
  marker_table_list[[i]] <- marker_table

  try({WriteXLS::WriteXLS(marker_table,
    figurePath(paste0(types[i],"/", "cluster_markers_", cell_names[i],".xlsx")),
    row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)})
}

length(marker_table_list)
# ==== mclapply for plotting cluster markers
for (i in seq_along(obj_list)) {
  print(names(obj_list[i]))
  mcPlotMarkers(obj_list[[i]], folder_prefix = paste0(types[i], "/"),
    marker_table_list[[i]], n_genes = 200)
}


