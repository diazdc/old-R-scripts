# Experimenting with different cell combinations

library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())
print(ls())

script_name <- paste0("ptime_subset_experiments_seurat3_v1.2")

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

central_subset <- which(cell_idents %in%
  c("central-cells") &
  data_sets %in% c("0min", "30min", "1hr", "1.5hr", "2hr", "3hr", "5hr"))

prog_subset <- which(cell_idents %in% c("HC-prog") &
  data_sets %in% c("3hr", "5hr")) 

HC_subset <- which(cell_idents %in% c("early-HCs", "mature-HCs") &
  data_sets %in% c("homeo"))

seurat_obj <- seurat_obj[,c(central_subset, prog_subset, HC_subset)]

to_remove <- readRDS(dataPath(
  "cells_to_remove_ptime_subset_central_seurat3_v1.0_.RDS"))
seurat_obj <- seurat_obj[, which(!colnames(seurat_obj) %in% to_remove)]

colnames(seurat_obj@meta.data)[10] <- "cell.type.old"


# ================================================================ Changing QC

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- subset(seurat_obj,
  subset = nFeature_RNA > 600 & nFeature_RNA < 1600)

test_obj <- seurat_obj
DefaultAssay(test_obj) <- "RNA"

dim(test_obj)
length(table(test_obj@assays$RNA@counts["slc1a3a",]))
length(table(test_obj@assays$RNA@data["slc1a3a",]))

seurat_obj <- NormalizeData(seurat_obj,
    normalization.method = "LogNormalize", scale.factor = 10000)

table(seurat_obj@meta.data$data.set)
summary(seurat_obj@meta.data$nFeature_RNA)
seurat_obj

png(figurePath("cell_QC.png"), width = 20, height = 10,
  units = "in", res = 300)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"),
  pt.size = 0.5, ncol = 2, group.by = "data.set")
dev.off()


# ================================================ Normalization/Dim Reduction


DefaultAssay(seurat_obj) <- "RNA"

if (FALSE) {
  seurat_obj <- SCTransform(seurat_obj)
} else {
  seurat_obj <- NormalizeData(seurat_obj,
    normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(
    seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
}

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


# ==== Plot by cell type and treatment
type_and_trt <- paste(seurat_obj@meta.data$data.set,
  seurat_obj@meta.data$cell.type.ident, sep = ".")

type_treat_levels <- c("0min.central-cells",
    "30min.central-cells", "1hr.central-cells", "1.5hr.central-cells",
    "2hr.central-cells", "3hr.central-cells", "5hr.central-cells",
    "1.5hr.HC-prog", "2hr.HC-prog", "3hr.HC-prog", "5hr.HC-prog",
    "10hr.HC-prog", "homeo.HC-prog", "3hr.early-HCs", "5hr.early-HCs",
    "10hr.early-HCs", "homeo.early-HCs", "3hr.mature-HCs", "5hr.mature-HCs",
    "10hr.mature-HCs", "homeo.mature-HCs")

seurat_obj@meta.data$cell.type.and.trt <- factor(type_and_trt,
  levels = type_treat_levels[type_treat_levels %in% type_and_trt],
  ordered = TRUE)

seurat_obj@meta.data[is.na(seurat_obj@meta.data$cell.type.and.trt),]
colors <- gg_color_hue(length(unique(seurat_obj@meta.data$cell.type.and.trt)))

p4 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.50,
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

if (FALSE) {
  testPCs(seurat_obj, from_to = 4:30, include_conditions = TRUE)
}

# ! Final Object
if(FALSE) {
  saveRDS(seurat_obj, dataPath(paste0("SeurObj_", script_name,"_.RDS")))
  seurat_obj <- readRDS(dataPath(paste0("SeurObj_", script_name,"_.RDS")))
}

if(FALSE) {
  source(paste0("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/",
    "Seurat3-final/ptime_subset_experiments_seurat3_v1.2_.R"), echo = TRUE)
}

stop("End Analysis")


# ========================================================== Additional Figures


# ==== mclapply for Find Markers
marker_table <- mcFindMarkers(seurat_obj,
  save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "clusters_markers_")

WriteXLS::WriteXLS(marker_table,
  figurePath(paste0("cluster_markers_", script_name,".xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)


# ==== mclapply for plotting cluster markers
mcPlotMarkers(seurat_obj, marker_table, n_genes = 200,
  folder_prefix = "cluster-markers-featureplot")

meta <- seurat_obj@meta.data
hr1_ind <- which(meta$data.set == "1hr")
meta$new_column <- ""
meta[hr1_ind,]$new_column <- paste(meta[hr1_ind,]$data.set,
  meta[hr1_ind,]$tree.ident, sep = ".")

seurat_obj@meta.data <- meta

png(figurePath("test.png"), width = 10, height = 10,
  units = "in", res = 300)
VlnPlot(seurat_obj, features = "slc1a3a",
  pt.size = 0.5, ncol = 1, group.by = "data.set")
dev.off()

table(seurat_obj@assays$RNA@counts["slc1a3a",])
table(seurat_obj@assays$RNA@data["slc1a3a",])

count_sum <- sum(seurat_obj@assays$RNA@counts["slc1a3a",])
count_vec <- seurat_obj@assays$RNA@counts["slc1a3a",]

nrml_vec <- vapply(count_vec, function(x) {
  scaled_dat <- ((x / count_sum) * 10000)
  }, double(1))

table(nrml_vec)