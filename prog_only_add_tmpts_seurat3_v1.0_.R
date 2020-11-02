# 1.0 Subsetted cells from all_LL_cells_regen_anchored_seurat3_v1.2_.R

library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("prog_only_add_tmpts_seurat3_v1.0")

dir.create(paste0("/n/projects/ddiaz/Analysis",
  "/Scripts/sb2191-regen/Seurat3-final/", script_name,"_figures/"),
  showWarnings = FALSE)

figurePath <- function(filename){paste0("/n/projects/ddiaz/Analysis",
  "/Scripts/sb2191-regen/Seurat3-final/", script_name, "_figures/", filename)}

dataPath <- function(filename){paste0("/n/projects/",
  "ddiaz/Analysis/Data/sb2191-regen/", filename)}

gene_info <- read.delim(paste0("/n/projects/ddiaz/Analysis/",
  "Data/gene-lists/Danio_Features_unique_Ens91_v2.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

trt_colors <- c("green3", "gold", "darkorange", "red", "magenta",
  "mediumpurple1", "lightseagreen", "deepskyblue", "blue")

ids <- c("homeo", "0min", "30min", "1hr", "1.5hr","2hr", "3hr", "5hr", "10hr")


# ============================================================ Reading in Data


seurat_obj <- readRDS(paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "SeurObj_anchored_cell_type_update_nmast_only_additional_timepoints",
  "_seurat3_v1.0_.RDS"))

ind <- which(seurat_obj@meta.data$cell.type.ident == "HC-prog")
seurat_obj <- seurat_obj[,ind]
seurat_obj@meta.data

DefaultAssay(seurat_obj) <- "RNA"


# ============================================================ Integration


seurat_obj_list <- SplitObject(seurat_obj, split.by = "data.set")
seurat_obj_list <- seurat_obj_list[ids]



dims = 1:15

SC_transform <- FALSE
if (SC_transform) {

  for (i in names(seurat_obj_list)) {
    seurat_obj_list[[i]] <- SCTransform(
      seurat_obj_list[[i]], verbose = FALSE)
  }
  seurat_obj_features <- SelectIntegrationFeatures(
    object.list = seurat_obj_list, nfeatures = 2000)
  
  seurat_obj_list <- PrepSCTIntegration(
    object.list = seurat_obj_list, anchor.features = seurat_obj_features)

  # reference = 1 is the homeostatic dataset in obj list
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
    anchor.features = seurat_obj_features, normalization.method = "SCT",
    dims = dims, reference = 1)

  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims,
    normalization.method = "SCT")
  
  } else {
  
  for (i in 1:length(seurat_obj_list)) {
    seurat_obj_list[[i]] <- NormalizeData(
      seurat_obj_list[[i]], verbose = FALSE)
    
    seurat_obj_list[[i]] <- FindVariableFeatures(
      seurat_obj_list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
    }
    
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
    dims = dims, reference = 1, k.filter = 84) # 1 is the homeostatic dataset in obj list

  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims)
  DefaultAssay(obj_integrated) <- "integrated"
}

obj_integrated@meta.data$data.set <- factor(
  obj_integrated@meta.data$data.set, ordered = TRUE, levels = ids)

if (FALSE) {
  saveRDS(obj_integrated, dataPath(
    paste0("SeurObj_before_clust", "_", script_name,"_.RDS")))

  obj_integrated <- readRDS(
    dataPath(paste0("SeurObj_before_clust", "_", script_name,"_.RDS")))
}


# =========================================================== UMAP/Clustering


if (!SC_transform) {
  obj_integrated <- ScaleData(obj_integrated, verbose = FALSE)
}

obj_integrated <- RunPCA(obj_integrated, npcs = 50, verbose = FALSE)

pdf(figurePath("PC_elbowPlot.pdf"))
ElbowPlot(obj_integrated)
dev.off()

if (FALSE) {
  DefaultAssay(obj_integrated) <- "integrated"
}

dims <- 1:7
res <- 0.2
obj_integrated <- FindNeighbors(obj_integrated, dims = dims, k.param = 20)
obj_integrated <- FindClusters(obj_integrated, resolution = res)
obj_integrated <- RunUMAP(obj_integrated, dims = dims,
  min.dist = 0.3, spread = 1)

obj_integrated <- BuildClusterTree(
  obj_integrated, reorder = TRUE, reorder.numeric = TRUE )


# ==== Plot by cluster
p1 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.60,
  label = TRUE, label.size = 4)
p1 <- cleanUMAP(p1)

png(figurePath(paste(
  "anchored_clusters_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p1)
dev.off()


# ==== Plot by treatment
p2 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.60,
  label = FALSE, label.size = 4, group.by = "data.set", cols = trt_colors)
p2 <- cleanUMAP(p2)

png(figurePath(paste("anchored_dataset_dims", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p2)
dev.off()


png(figurePath(paste0("UMI_count_by_dataset.png")),
  width = 12, height = 10, units = "in", res = 300)
VlnPlot(obj_integrated,"nCount_RNA",
  group.by= "data.set", pt.size = 0.3)
dev.off()


fplot <- FeaturePlot(obj_integrated,"nCount_RNA",pt.size = 0.6)
fplot <- cleanUMAP(fplot)

png(figurePath(paste0("UMI_count_.png")),
  width = 12, height = 10, units = "in", res = 300)
print(fplot)
dev.off()

stop("End Analysis")

if(FALSE) {
  source(paste0("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/",
    "Seurat3-final/prog_only_add_tmpts_seurat3_v1.0_.R"),
    echo = TRUE)
}

if (FALSE) {
  testPCs(obj_integrated, from_to = 16:25)
}


# ========================================================== Additional Figures


DefaultAssay(obj_integrated) <- "RNA"

# ==== mclapply for Find Markers
marker_table <- mcFindMarkers(obj_integrated,
  save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "res0.2_")

WriteXLS::WriteXLS(marker_table,
  figurePath(paste0("cluster_markers_res0.2", script_name,".xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

# ==== mclapply for plotting cluster markers
mcPlotMarkers(obj_integrated, marker_table, n_genes = 200,
  folder_prefix = "cluster-markers-res0.2")



