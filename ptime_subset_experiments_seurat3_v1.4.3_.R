# Experimenting with different cell combinations

{
library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())
print(ls())

script_name <- paste0("ptime_subset_experiments_seurat3_v1.4.3")

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
}


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
  data_sets %in% c(
    "homeo", "3hr", "5hr"))

prog_subset <- which(cell_idents %in% c("HC-prog") &
  data_sets %in% c("homeo", "1hr", "1.5hr", "2hr", "3hr", "5hr", "10hr"))

early_HC_subset <- which(cell_idents %in% c("early-HCs") &
  data_sets %in% c("1.5hr", "2hr", "3hr", "5hr", "10hr"))

HC_subset <- which(cell_idents %in% c("mature-HCs") &
  data_sets %in% c("homeo", "2hr", "3hr", "5hr", "10hr"))

seurat_obj <- seurat_obj[,
  c(central_subset, prog_subset, early_HC_subset, HC_subset)]

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

type_treat_levels <- c("homeo.central-cells", "0min.central-cells",
  "30min.central-cells", "1hr.central-cells", "1.5hr.central-cells",
  "2hr.central-cells", "3hr.central-cells", "5hr.central-cells",
  "30min.HC-prog", "1hr.HC-prog", "1.5hr.HC-prog", "2hr.HC-prog", "3hr.HC-prog",
  "5hr.HC-prog", "10hr.HC-prog", "homeo.HC-prog", "1.5hr.early-HCs",
  "2hr.early-HCs", "3hr.early-HCs", "5hr.early-HCs", "10hr.early-HCs",
  "homeo.early-HCs", "1.5hr.mature-HCs", "2hr.mature-HCs", "3hr.mature-HCs",
  "5hr.mature-HCs", "10hr.mature-HCs", "homeo.mature-HCs")

seurat_obj@meta.data$cell.type.and.trt <- factor(type_and_trt,
  levels = type_treat_levels[type_treat_levels %in% type_and_trt],
  ordered = TRUE)

seurat_obj@meta.data[is.na(seurat_obj@meta.data$cell.type.and.trt),]
colors <- gg_color_hue(length(unique(seurat_obj@meta.data$cell.type.and.trt)))
n_levels_cell_trt <- length(levels(seurat_obj@meta.data$cell.type.and.trt))

# Stagger colors for more contrast
colors <- colors[c(seq(1, n_levels_cell_trt, by = 2),
  seq(2,n_levels_cell_trt, by = 2))]

# Other staggered color option
colors2 <- gg_color_hue(
  length(unique(seurat_obj@meta.data$cell.type.and.trt)) / 3)
colors2 <- rep(colors2, 3)

# All grey cells except chosen cell type
ind <- grep("central", unique(seurat_obj@meta.data$cell.type.and.trt))
cols_short <- gg_color_hue(length(ind))
greys <- rep("#E4E4E4", length(unique(seurat_obj@meta.data$cell.type.and.trt)))
greys[ind] <- cols_short
colors3 <- greys

p4 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.50,
  group.by = "cell.type.and.trt", label = FALSE, label.size = 4,
  cols = colors3) + theme(legend.position="bottom")
p4 <- cleanUMAP(p4)

png(figurePath(paste(
  "cell_type_and_trt_central_", "dim", max(dims),".png", sep = "_")),
  width = 11, height = 11, units = "in", res = 300)
print(p4)
dev.off()

p5 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.50,
  group.by = "cell.type.and.trt", label = FALSE, label.size = 4,
  cols = colors) + theme(legend.position="bottom")
p5 <- cleanUMAP(p5)

png(figurePath(paste(
  "cell_type_and_trt_", "dim", max(dims),".png", sep = "_")),
  width = 11, height = 11, units = "in", res = 300)
print(p5)
dev.off()

table(seurat_obj@meta.data$cell.type.and.trt)

# UMI counts for each dataset
png(figurePath(paste0("UMI_count_by_dataset.png")),
  width = 12, height = 12, units = "in", res = 300)
VlnPlot(seurat_obj,"nCount_RNA",
  group.by= "data.set", pt.size = 0.3)
dev.off()

if (FALSE) {
  testPCs(seurat_obj, from_to = 4:15, include_conditions = TRUE)
}

# ! Final Object
if(FALSE) {
  saveRDS(seurat_obj, dataPath(paste0("SeurObj_", script_name,"_.RDS")))
  seurat_obj <- readRDS(dataPath(paste0("SeurObj_", script_name,"_.RDS")))
}

if(FALSE) {
  source(paste0("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/",
    "Seurat3-final/ptime_subset_experiments_seurat3_v1.4.3_.R"), echo = TRUE)
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


# ==== LDA topic analysis
library(CountClust)

if (FALSE) {
# Prob. LDA clustering (TAKES FOREVER and uses a TON of RAM) 10 (Rsession2),
# 15 (Rsession 3), 20 (Rsession 4), 25 (session 5), 30 (Rsession6)
# Session 4 and 5 are on cedar. all others on my linux box

obj_counts <- as.matrix(seurat_obj@assays$RNA@data)
obj_meta <- seurat_obj@meta.data
gene_names <- rownames(obj_counts)

start_time <- Sys.time()
obj_FitGoM <- FitGoM(t(obj_counts), K = 15)
end_time <- Sys.time()
end_time - start_time

if(FALSE) {
  saveRDS(obj_FitGoM, dataPath(paste0("FitGoM_obj_topics_",
    obj_FitGoM$fit$K, ".RDS")))

  n_topics <- 15
  obj_FitGoM <- readRDS(dataPath(paste0("FitGoM_obj_topics_",
    n_topics, ".RDS")))
}


Idents(seurat_obj) <- seurat_obj@meta.data$cell.type.and.trt
omega <- data.frame(obj_FitGoM$fit$omega)
annotation <- data.frame(sample_id = rownames(omega),
  tissue_label = paste0("", seurat_obj@active.ident))

colnames(omega) <- paste0("topic", 1:obj_FitGoM$fit$K)
rownames(omega) <- annotation$sample_id;
annotation$tissue_label <- as.factor(annotation$tissue_label)

GoM_plot <- StructureGGplot(omega = omega,
  annotation = annotation, palette = gg_color_hue(obj_FitGoM$fit$K),
  yaxis_label = "Cells", order_sample = TRUE,
  axis_tick = list(axis_ticks_length = .1, axis_ticks_lwd_y = .1,
  axis_ticks_lwd_x = .1, axis_label_size = 7, axis_label_face = "bold"))

png(figurePath(paste0("GoM_topic_",obj_FitGoM$fit$K,"_.png")),
  width = 8, height = 10, res = 300, units = "in")
print(GoM_plot)
dev.off()

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, omega)
seurat_obj <- AddMetaData(seurat_obj, omega)

to_plot <- unlist(lapply(1:obj_FitGoM$fit$K,
  function(i){
    paste0("topic",i)
}))

topic_umap <- FeaturePlot(object = seurat_obj,
  features = to_plot, combine = FALSE)

for (i in 1:length(topic_umap)) {
  topic_umap[[i]] <- cleanUMAP(topic_umap[[i]])
  topic_umap[[i]] <- topic_umap[[i]] + ggtitle(paste0("topic", i))
}

png(figurePath(paste0("umap_topic_",obj_FitGoM$fit$K,"_.png")),
  width = 25, height = 40, res = 300, units = "in")
print(cowplot::plot_grid(plotlist = topic_umap, ncol = 3))
dev.off()

# Generate table with genes
theta_mat <- obj_FitGoM$fit$theta
top_features <- ExtractTopFeatures(theta_mat, top_features = 100,
  method = "poisson", options = "min")

gene_list <- do.call(rbind, lapply(1:dim(top_features$indices)[1],
  function(x) gene_names[top_features$indices[x,]]))

topic_gene_table <- as.data.frame(t(gene_list))
colnames(topic_gene_table) <- paste0(
  rep("topic", obj_FitGoM$fit$K), seq(1, obj_FitGoM$fit$K))

WriteXLS::WriteXLS(topic_gene_table,
  figurePath(paste0("topic_", obj_FitGoM$fit$K,"_gene_table.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

}



# ==== Transpose data for gene module detection
{
n_genes <- 2000
aux_obj <- FindVariableFeatures(
  seurat_obj, selection.method = "vst", nfeatures = n_genes)
genes_to_keep <- aux_obj@assays$RNA@var.features

cnt_mtx <- t(as.matrix(seurat_obj@assays$RNA@counts))
gene_obj <- CreateSeuratObject(cnt_mtx)
gene_obj <- gene_obj[,which(colnames(gene_obj) %in% genes_to_keep)]
gene_obj <- NormalizeData(gene_obj)

# nfeatures are now cells, do not remove any.
gene_obj <- FindVariableFeatures(
  gene_obj, selection.method = "vst", nfeatures = nrow(gene_obj))
gene_obj <- ScaleData(gene_obj, features = rownames(gene_obj))

gene_obj <- RunPCA(gene_obj, features = VariableFeatures(object = gene_obj))
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

pdf(figurePath("gene_PC_elbowPlot.pdf"))
ElbowPlot(seurat_obj)
dev.off()

res <- 0.4
dims <- 1:11
gene_obj <- FindNeighbors(gene_obj, dims = dims)
gene_obj <- FindClusters(gene_obj, resolution = res)
gene_obj <- RunUMAP(gene_obj, dims = dims)

gene_obj <- BuildClusterTree(
  gene_obj, reorder = TRUE, reorder.numeric = TRUE )

# Plot by cluster
p5 <- DimPlot(gene_obj, reduction = "umap", pt.size = 0.40,
  label = TRUE, label.size = 4)
p5 <- cleanUMAP(p5)

testPCs(gene_obj, from_to = 16:25, include_conditions = FALSE)

png(figurePath(paste(
  "gene_modules_", res, "dim", max(dims), n_genes,"_genes_.png", sep = "_")),
  width = 10, height = 8, units = "in", res = 300)
print(p5)
dev.off()

module_table <- data.frame(Gene.name.uniq = rownames(gene_obj@meta.data),
  module = gene_obj@meta.data$tree.ident)

module_table <- module_table[order(module_table$module),]
WriteXLS::WriteXLS(module_table,
  figurePath(paste0("gene_module_table_", n_genes,"_genes_.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)
}