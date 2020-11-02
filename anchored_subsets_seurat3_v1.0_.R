library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("anchored_subsets_seurat3_v1.0")

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

ids <- c("homeo", "0min", "30min", "1hr", "1.5hr","2hr", "3hr", "5hr", "10hr")


# ============================================================ Reading in Data


seurat_obj <- readRDS(paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "SeurObj_anchored_cell_type_update_nmast_only_additional_timepoints",
  "_seurat3_v1.0_.RDS"))

seurat_obj@meta.data$cell.type.ident <- plyr::revalue(
  seurat_obj@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))

DefaultAssay(seurat_obj) <- "RNA"

# Select cell types.
meta <- seurat_obj@meta.data
types <- levels(meta$cell.type.ident)[c(4:6,8)]
cell_names <- types
obj_unint_list <- list()[1:4]

for (i in 1:4){
  obj_unint_list[[i]] <- seurat_obj[,meta$cell.type.ident == types[i]]
}

obj_unint_list <- setNames(obj_unint_list, types[1:4])


# ============================================================== Integration


seurat_obj_list <- list()[1:4]
for (i in 1:4) {
  seurat_obj_list[[i]] <- SplitObject(obj_unint_list[[i]],
    split.by = "data.set")
}

seurat_obj_list <- setNames(seurat_obj_list , types)

# print cell counts
for (i in 1:4) {
  print(names(seurat_obj_list[i]))
  print(vapply(seq_along(ids), function(k) {
      ncol(seurat_obj_list[[i]][[k]])
    }, integer(1)))
}

seurat_obj_list <- setNames(seurat_obj_list, types[1:4])
obj_list <- list()[1:4]
for (i in 1:4) {
  for (j in 1:length(seurat_obj_list[[i]])) {
    seurat_obj_list[[i]][[j]] <- NormalizeData(
      seurat_obj_list[[i]][[j]], verbose = FALSE)
    
    seurat_obj_list[[i]][[j]] <- FindVariableFeatures(
      seurat_obj_list[[i]][[j]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
    seurat_obj_list[[i]][[j]]
  }
  print(i)
  cnts <- vapply(seq_along(ids), function(k) {
      ncol(seurat_obj_list[[i]][[k]])
    }, integer(1))

  print(names(seurat_obj_list[i]))
  print(names(seurat_obj_list[[i]]))
  print((cnts))
  seq_along(min(cnts - 1))

  obj_anchors <- FindIntegrationAnchors(
    object.list = seurat_obj_list[[i]], dims = 1:15,
    k.filter = min(cnts), reference = 1)

  obj_list[[i]] <- IntegrateData(anchorset = obj_anchors,
    dims = 1:15)
  DefaultAssay(obj_list[[i]]) <- "integrated"

  obj_list[[i]] <- ScaleData(obj_list[[i]], verbose = FALSE)

  obj_list[[i]]@meta.data$data.set <- factor(
    obj_list[[i]]@meta.data$data.set, ordered = TRUE, levels = ids)
}


# ============================================================== Dim reduction


for (i in 1:4){
  dir.create(figurePath(types[i]), showWarnings = FALSE)
}

for (i in 1:4) {
  obj_list[[i]] <- RunPCA(obj_list[[i]],
    features = VariableFeatures(obj_list[[i]]))

  pc_dims <- 1:50
  if (TRUE) {
      pdf(figurePath(paste0(types[i],"/",
        cell_names[i], "_PC_elbowPlot.pdf")))
      ElbowPlot(obj_list[[i]], ndims = 50)
  }
}


# ======================================================== Clustering and UMAP


dims <- 1:8
res = 0.1

for (i in 1:4) {
  obj_list[[i]] <- FindNeighbors(obj_list[[i]], dims = dims, k.param = 20)
  obj_list[[i]] <- RunUMAP(obj_list[[i]],
    dims = dims, min.dist = 0.3, spread = 1)
}

for (i in 1:4) {
  obj_list[[i]] <- FindClusters(obj_list[[i]], resolution = res)
  obj_list[[i]] <- BuildClusterTree(obj_list[[i]],
    reorder = TRUE, reorder.numeric = TRUE)
}


# ==== Plot by cluster
for (i in 1:4) {
  p1 <- DimPlot(obj_list[[i]], reduction = "umap", pt.size = 0.40,
    label = TRUE, label.size = 4)
  p1 <- cleanUMAP(p1)

  png(figurePath(paste0(types[i],"/", cell_names[i], "_UMAP_clusters_",
    max(dims), "_", res, ".png")), width = 10, height = 8, units = "in",
    res = 300)
  print(p1)
  dev.off()
}


# ==== Plot by treatment
for (i in 1:4) {
  p2 <- DimPlot(obj_list[[i]], reduction = "umap", pt.size = 0.40,
    label = FALSE, label.size = 4, group.by = "data.set", cols = trt_colors)
  p2 <- cleanUMAP(p2)

  png(figurePath(paste0(types[i],"/", cell_names[i],
    "_UMAP_dataset_", max(dims),".png")),
    width = 10, height = 8, units = "in", res = 300)
  print(p2)
  dev.off()
}


# ============================================================= Save/load data


if(TRUE) {
  saveRDS(obj_list,
    dataPath(paste0("CombinedObj_", script_name,"_.RDS")))
  obj_list <- readRDS(
    dataPath(paste0("CombinedObj_", script_name,"_.RDS")))
}

for (i in 1:4) {
  str(obj_list[[i]]@assays$RNA@counts)
}

table(seurat_obj@meta.data$data.set)
stop("break section 2")


# ========================================================= Additional figures


for (i in 1:4) {
  DefaultAssay(obj_list[[i]]) <- "RNA"
}

# ==== mclapply for Find Markers
marker_table_list <- list()[1:4]
for (i in 1:4) {
  marker_table <- mcFindMarkers(obj_list[[i]],
    save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "")
  marker_table_list[[i]] <- marker_table

  WriteXLS::WriteXLS(marker_table,
    figurePath(paste0(types[i],"/", "cluster_markers_", script_name,".xlsx")),
    row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)
}


# ==== mclapply for plotting cluster markers
for (i in 1:4) {
mcPlotMarkers(obj_list[[i]], marker_table_list[[i]], n_genes = 200,
  folder_prefix = paste0(types[i],"/"))
}


# ==== mclapply for plotting cluster markers
for (i in 1:4) {
mcPlotMarkers(obj_list[[i]], marker_table_list[[i]], n_genes = 200,
  folder_prefix = paste0(types[i],"/"))
}


# ==== mclapply for plotting cluster markers
for (i in 1:4) {
  p3 <- FeaturePlot(obj_list[[i]], features = "nCount_RNA", pt.size = 0.40)
  p3 <- cleanUMAP(p3)

  png(figurePath(paste0(types[i],"/", cell_names[i],
    "_UMAP_nUMI_", max(dims),".png")),
    width = 10, height = 8, units = "in", res = 300)
  print(p3)
  dev.off()
}
