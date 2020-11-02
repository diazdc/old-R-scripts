# Removes homeo progentiors

library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("ptime_subset_experiments_seurat3_v1.1")

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


seurat_obj <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
    "SeurObj_anchored_cell_type_update_nmast_",
    "only_additional_timepoints_seurat3_v1.0_.RDS"))

cell_idents <- seurat_obj@meta.data$cell.type.ident
data_sets <- seurat_obj@meta.data$data.set
levels(cell_idents)

central_subset <- which(cell_idents %in% c("central-cells", "HC-prog") &
  data_sets %in% c("3hr", "5hr"))

prog_subset <- which(cell_idents %in% c("HC-prog", "early-HCs", "mature-HCs") &
  data_sets %in% c("homeo"))

seurat_obj <- seurat_obj[,c(central_subset, prog_subset)]
colnames(seurat_obj@meta.data)[10] <- "cell.type.old"


# ================================================ Normalization/Dim Reduction


DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- NormalizeData(seurat_obj,
  normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- FindVariableFeatures(
  seurat_obj, selection.method = "vst", nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

pdf(figurePath("PC_elbowPlot.pdf"))
ElbowPlot(seurat_obj)
dev.off()

dims <- 1:9
res <- 0.6
seurat_obj <- FindNeighbors(seurat_obj, dims = dims, k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = res)
seurat_obj <- RunUMAP(seurat_obj, dims = dims,
  min.dist = 0.3, spread = 1)

seurat_obj <- BuildClusterTree(
  seurat_obj, reorder = TRUE, reorder.numeric = TRUE )


# ==== Plot by cluster
p1 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.40,
  label = TRUE, label.size = 4)
p1 <- cleanUMAP(p1)

png(figurePath(paste(
  "clusters_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p1)
dev.off()


# ==== Plot by treatment
p2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.40,
  label = FALSE, label.size = 4, group.by = "data.set", cols = trt_colors)
p2 <- cleanUMAP(p2)

png(figurePath(paste("dataset_dims", max(dims),".png", sep = "_")),
  width = 11, height = 9, units = "in", res = 300)
print(p2)
dev.off()


# ==== Plot by cell type
colors <- gg_color_hue(length(levels(seurat_obj@meta.data$cell.type.ident)))

p3 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.40,
  group.by = "cell.type.ident", label = TRUE, label.size = 4, cols = colors)
p3 <- cleanUMAP(p3)

png(figurePath(paste(
  "cell_type", "dim", max(dims),".png", sep = "_")),
  width = 11, height = 9, units = "in", res = 300)
print(p3)
dev.off()


# ==== Plot by cell type
seurat_obj[["cell.type.and.trt"]] <- factor(paste(
  seurat_obj@meta.data$data.set, seurat_obj@meta.data$cell.type.ident,
  sep = "."))

colors <- gg_color_hue(length(levels(seurat_obj@meta.data$cell.type.and.trt)))

p4 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.60,
  group.by = "cell.type.and.trt", label = FALSE, label.size = 4, cols = colors)
p4 <- cleanUMAP(p4)

png(figurePath(paste(
  "cell_type_and_trt", "dim", max(dims),".png", sep = "_")),
  width = 11, height = 9, units = "in", res = 300)
print(p4)
dev.off()


# UMI counts for each dataset
png(figurePath(paste0("UMI_count_by_dataset.png")),
  width = 12, height = 10, units = "in", res = 300)
VlnPlot(seurat_obj,"nCount_RNA",
  group.by= "data.set", pt.size = 0.3)
dev.off()

stop("End Analysis")

if(FALSE) {
  source(paste0("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/",
    "Seurat3-final/ptime_subset_experiments_seurat3_v1.0_.R"),
    echo = TRUE)
}

if (FALSE) {
  testPCs(seurat_obj, from_to = 5:20, include_conditions = TRUE)
}

png(figurePath(paste0("test_plot.png")),
  width = 12, height = 10, units = "in", res = 300)
print(DotPlot(seurat_obj, features = c("atoh1a", "dld", "dla")))
dev.off()


# ==== median gene count per cluster
vln_list <- list()[1:9]
vln_list <- lapply(seq_along(ids), function(i) {
  ind <- which(obj_integrated@meta.data$data.set == ids[i])
  VlnPlot(obj_integrated[,ind],"nCount_RNA",
    group.by= "cell.type.ident", pt.size = 0.3) + ggtitle(ids[i])
  })

png(figurePath(paste0("UMI_count_by_cell_type.png")),
  width = 25, height = 15, units = "in", res = 300)
print(cowplot::plot_grid(plotlist = vln_list), ncol = 3)
dev.off()


# ! Final Object
if(TRUE) {
  saveRDS(seurat_obj, dataPath(paste0("SeurObj2_", script_name,"_.RDS")))
  seurat_obj <- readRDS(dataPath(paste0("SeurObj2_", script_name,"_.RDS")))
}

