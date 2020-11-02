# 1.0 Subsetted cells from all_LL_cells_regen_anchored_seurat3_v1.2_.R

{
library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("additional_timepoints_seurat3_v1.0")

dir.create(paste0("/n/projects/ddiaz/Analysis",
  "/Scripts/sb2191-regen/Seurat3-final/", script_name,"_figures/"))

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

# trt_colors <- gg_color_hue(9)

trt_colors <- c("green3", "gold", "darkorange", "red", "magenta",
  "mediumpurple1", "lightseagreen", "deepskyblue", "blue")

ids <- c("homeo", "0min", "30min", "1hr", "1.5hr","2hr", "3hr", "5hr", "10hr")
}

# ============================================================ Reading in Data


paths <- Sys.glob(paste0("/n/projects/ddiaz/Analysis/Data/MOLNG*",
  "/raw_feature_bc_matrix*/zipped/"))

# ordering to timepoint, omitting 10hr (2588) small dataset
paths <- paths[c(1,8,7,6,9,10,2:5)]

raw_mats <- list()[9]
for (i in 1:9){
  raw_mats[[i]] <- Read10X(paths[i])
}

all_objs <- list()[9]

for (i in 1:9){
  all_objs[[i]] <- CreateSeuratObject(raw_mats[[i]], project = ids[i],
    min.cells = 5, min.features = 300)
}

for (i in 1:9) {
  print(ids[i])
  print("genes cells")
  print(dim(all_objs[[i]]))
  cat("\n")
}

seurat_obj <- merge(all_objs[[1]], y = all_objs[2:9], add.cell.ids = ids)
nrow(seurat_obj@meta.data)

seurat_obj@meta.data$data.set <- factor(sub(
  "_(.+)", "",rownames(seurat_obj@meta.data)), levels = ids, ordered = TRUE)

table(seurat_obj@meta.data$data.set)


# =============================================================== Downsampling


DefaultAssay(seurat_obj) <- "RNA"
table(seurat_obj@meta.data$data.set)

# Downsampling the 1hr data set
avg <- ceiling(mean(table(seurat_obj@meta.data$data.set)[2:9]))
avg_ind_1hr <- sample(grep(paste0("^", "1hr"), colnames(seurat_obj)), avg)

other_ind <- grep("^((?!1hr).)*$", colnames(seurat_obj), perl = TRUE)
seurat_obj_sub <- seurat_obj[ ,c(other_ind, avg_ind_1hr)]

# Check cell numbers
table(seurat_obj_sub@meta.data$data.set)
sum(table(seurat_obj_sub@meta.data$data.set))


# Section 1 Anchoring Cells
# ================================================================ Integration


seurat_obj_list <- SplitObject(seurat_obj_sub, split.by = "data.set")
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
    dims = dims, reference = 1) # 1 is the homeostatic dataset in obj list

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

if (FALSE){
  DefaultAssay(obj_integrated) <- "integrated"
}

dims <- 1:34
res <- 1.8
obj_integrated <- FindNeighbors(obj_integrated, dims = dims, k.param = 20)
obj_integrated <- FindClusters(obj_integrated, resolution = res)
obj_integrated <- RunUMAP(obj_integrated, dims = dims,
  min.dist = 0.3, spread = 1)

obj_integrated <- BuildClusterTree(
  obj_integrated, reorder = TRUE, reorder.numeric = TRUE )


# ==== Plot by cluster
p1 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.20,
  label = TRUE, label.size = 4)
p1 <- cleanUMAP(p1)

png(figurePath(paste(
  "anchored_clusters_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p1)
dev.off()


# ==== Plot by treatment
p2 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.20,
  label = FALSE, label.size = 4, group.by = "data.set", cols = trt_colors)
p2 <- cleanUMAP(p2)

png(figurePath(paste("anchored_dataset_dims", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p2)
dev.off()

if (FALSE){
  testPCs(obj_integrated, from_to = 10:50)
}


# ============================================================= Save/load data


if(TRUE) {
  saveRDS(obj_integrated,
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
  print("Object saved")

  obj_integrated <- readRDS(
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
}

stop("Analysis complete")

table(obj_integrated@meta.data$data.set)
table(obj_integrated@meta.data$tree.ident)

# Adding in cell type names
meta <- obj_integrated@meta.data
colnames(meta)

cells <- list("mature-HCs" = 1, "early-HCs" = c(19,24),  "HC-prog" = 17,
  "central-cells" = c(6:9), "DV-cells" = 4, "AP-cells" = c(2,3),
  "amp-SCs" = c(5,13,18), "mantle-cells" = c(14,15), "inm" = c(26:33),
  "col1a1b-pos" = c(21,34), "c1qtnf5-pos" = 20, "clec14a-pos" = 25,
  "laptm5-pos" = 22, "foxi3b-pos" = 23, "cldn1-pos" = 11, "ela2l-pos" = 16,
  "hbbe1.1-pos" = 12, "agxtb-pos" = 10)

meta$cell.type.ident <- factor(rep("", nrow(meta)),
  levels = names(cells), ordered = TRUE)

for (i in 1:length(cells)) {
  meta$cell.type.ident[meta$tree.ident %in% cells[[i]]] <- names(cells)[i]
}

gen_type <- list()[seq_along(cells)]
gen_type[seq_along(cells)] <- names(cells)
names(gen_type) <- c(rep("neuromast", 8), rep("non-specific", 8))

meta$general.ident <- factor(rep("", nrow(meta)),
  levels = unique(names(gen_type)), ordered = TRUE)

for (i in 1:length(cells)) {
  meta$general.ident[
    meta$cell.type.ident %in% gen_type[[i]]] <- names(gen_type)[i]
}

obj_integrated@meta.data <- meta
Idents(obj_integrated) <- obj_integrated@meta.data$cell.type.ident

p3 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.20,
  label = TRUE, label.size = 4)
p3 <- cleanUMAP(p3)

png(figurePath(paste(
  "anchored_cell_type_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p3)
dev.off()

# ! Final Object 
if(TRUE) {
  saveRDS(obj_integrated, dataPath(paste0(
    "SeurObj_anchored_cell_type_update", script_name,"_.RDS")))
  obj_integrated <- readRDS(dataPath(paste0(
    "SeurObj_anchored_cell_type_update", script_name,"_.RDS")))
}


# ========================================================= Additional figures
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
DefaultAssay(obj_integrated) <- "RNA"


# ==== mclapply for Find Markers
marker_table <- mcFindMarkers(obj_integrated,
  save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "cell_type_")

WriteXLS::WriteXLS(marker_table,
  figurePath(paste0("cell_type_markers_", script_name,".xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

marker_table <- readRDS(dataPath(paste0(
  "Clusters_anchored_DF", script_name,"_.RDS")))


# ==== mclapply for plotting
mcPlotMarkers(obj_integrated, marker_table, n_genes = 400,
  folder_prefix = "cell-type")


# ==== treatments by diff exp
levels(Idents(obj_integrated))
single_group <- diffConditionMrkrs(obj_integrated, n_cores = 1,
  group_clusters = c("AP-cells"), file_prefix = "AP")

WriteXLS::WriteXLS(single_group,
  figurePath(paste0("AP_combined_vs_other_all_timepoints.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

mantle_amp_AP <- diffConditionMrkrs(obj_integrated, n_cores = 1,
  group_clusters = c("mantle-cells", "amp-SCs", "AP-cells"))

WriteXLS::WriteXLS(mantle_amp_AP,
  figurePath(paste0("mantle_amp_AP_combined_vs_other_all_timepoints.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

pbulk <- diffConditionMrkrs(obj_integrated, n_cores = 1,
  cell_group_name = "pbulk", group_clusters = levels(Idents(obj_integrated)))

WriteXLS::WriteXLS(pbulk,
  figurePath(paste0("pbulk_all_timepoints.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)


# ==== plot treatment by diff exp
central_markers <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
    "central_all_markers_additional_timepoints_seurat3_v1.0_.RDS"))

central_markers <- central_markers[order(central_markers$p_val_adj,
  central_markers$cell.type.and.trt),]

diffConditionPlots(obj_integrated, input_file = central_markers,
  n_cores = 10, n_genes = 400, folder_prefix = "central-all-timepoints")

central_prog <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
    "central_and_prog_all_markers_additional_timepoints_seurat3_v1.0_.RDS"))

central_prog <- central_prog[order(central_prog$p_val_adj,
  central_prog$cell.type.and.trt),]

diffConditionPlots(obj_integrated, central_prog, n_cores = 10, n_genes = 400)

pbulk <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
    "pbulk_all_markers_additional_timepoints_seurat3_v1.0_.RDS"))

pbulk <- pbulk[order(pbulk$p_val_adj, pbulk$cell.type.and.trt),]

diffConditionPlots(obj_integrated, pbulk, all_idents = TRUE,
  folder_prefix = "pbulk", n_cores = 10, n_genes = 400)


# ==== Correlated genes
cor_list <- list()

sox4a_cor_list <- lapply(1:9, function(i) {
    cor_list[[i]] <- corrleatedGenes(obj_integrated, "sox4a*1",
    condition = ids[i], n_results = 100, verbose = FALSE)
    colnames(cor_list[[i]]) <- c(
      paste0(ids[[i]], "_gene"), paste0(ids[[i]], "_coefficient"))
    return(cor_list[[i]])
  }
)

sox4a_cor_df <- dplyr::bind_cols(sox4a_cor_list)
write.table(sox4a_cor_df, figurePath(paste0("sox4a_correlates.tsv")),
  row.names = FALSE, col.names = TRUE, sep = "\t")


# ==== Subsettable heatmaps
if (TRUE) {
  DefaultAssay(obj_integrated) <- "RNA"
  assay <- "assay_RNA_"
}

obj_integrated <- ScaleData(obj_integrated)
all_markers <- readRDS(dataPath(paste0("AllMarkers_", script_name,"_.RDS")))

# Subset by identity
meta <- obj_integrated@meta.data

treatment <- "central-cells"
BCs_to_subset <- rownames(meta[meta$cell.type.ident == treatment,])

# Subset object
obj_subset <- obj_integrated[,BCs_to_subset]

# Order genes by logFC
all_markers <- all_markers[order(all_markers$cell.type.and.trt,
  all_markers$avg_logFC, decreasing = TRUE),]
head(all_markers)

# Grab genes
# Display all possible comparisions
unique(all_markers$cell.type.and.trt)
cell_type_trt <- "central-cells_3hr"

to_plot <- all_markers$Gene.name.uniq[
  all_markers$cell.type.and.trt == cell_type_trt]
to_plot <- to_plot[1:100] # Number of genes to plot

hmap <- DoHeatmap(obj_subset, features = to_plot, slot = "scale.data",
  group.by = "data.set") + ggplot2::ggtitle(cell_type_trt) + 
  ggplot2::theme(plot.title = ggplot2::element_text(
    hjust = 0.5, size = 16, face = "bold")) 

hmap <- hmap + viridis::scale_fill_viridis()

png(figurePath(paste0("hmap_test.png")), width = 14, height = 14,
  units = "in", res = 200)
print(hmap)
dev.off()


# ==== Plot common features
DefaultAssay(obj_integrated) <- "RNA"

common_features <- scan(paste0("/n/projects/ddiaz/Analysis/",
  "Data/gene-lists/common_features.txt"), what = "character")

e <- FeaturePlot(obj_integrated, common_features,
  reduction = "umap", pt.size = 0.25, combine = FALSE)

for (i in 1:length(e)) {
  e[[i]] <- e[[i]] + NoLegend() + NoAxes()
}

png(figurePath("common_features.png"), width = 40,
  height = 80, units = "in", res = 200)
print(cowplot::plot_grid(plotlist = e, ncol = 4))
dev.off()

# source("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/Seurat3-final/additional_timepoints_seurat3_v1.0_.R", echo = TRUE)
