# 1.0 Subsetted cells from all_LL_cells_regen_anchored_seurat3_v1.2_.R

{
library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("nmast_skin_additional_timepoints_seurat3_v1.0")

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


obj_integrated <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
     "SeurObj_anchored_cell_type_updateadditional_",
     "timepoints_seurat3_v1.0_.RDS"))

to_subset <- obj_integrated@meta.data$cell.type.ident %in% 
  levels(Idents(obj_integrated))[c(1:8,10:12)]

obj_integrated <- obj_integrated[,to_subset]
colnames(obj_integrated@meta.data)[10] <- "cell.type.old"

table(obj_integrated@meta.data$nCount_RNA < 15000)
obj_integrated <- obj_integrated[,obj_integrated@meta.data$nCount_RNA < 15000]


# ============================================================ Integration


seurat_obj_list <- SplitObject(obj_integrated, split.by = "data.set")
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

dims <- 1:10
res <- 0.6
obj_integrated <- FindNeighbors(obj_integrated, dims = dims, k.param = 20)
obj_integrated <- FindClusters(obj_integrated, resolution = res)
obj_integrated <- RunUMAP(obj_integrated, dims = dims,
  min.dist = 0.3, spread = 1)

obj_integrated <- BuildClusterTree(
  obj_integrated, reorder = TRUE, reorder.numeric = TRUE )


# ==== Plot by cluster
p1 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.30,
  label = TRUE, label.size = 4)
p1 <- cleanUMAP(p1)

png(figurePath(paste(
  "anchored_clusters_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p1)
dev.off()


# ==== Plot by treatment
p2 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.30,
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

stop("End Analysis")

if(FALSE) {
  source(paste0("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/",
    "Seurat3-final/nmast_skin_additional_timepoints_seurat3_v1.0_.R"),
    echo = TRUE)
}

if (FALSE){
  testPCs(obj_integrated, from_to = 5:9)
}


# ============================================================= Save/load data


if(FALSE) {
  saveRDS(obj_integrated,
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
  print("Object saved")

  obj_integrated <- readRDS(
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
}

table(obj_integrated@meta.data[
  obj_integrated@meta.data$data.set == "homeo", "tree.ident"])
table(obj_integrated@meta.data$tree.ident)

# check res first. using res 0.6 for cell classification
levels(Idents(obj_integrated))
colnames(obj_integrated@meta.data)
Idents(obj_integrated) <- obj_integrated@meta.data$integrated_snn_res.0.6

# Adding in cell type names
meta <- obj_integrated@meta.data
colnames(meta)

cells <- list("mature-HCs" = 1, "early-HCs" = 2,  "HC-prog" = 8,
  "central-cells" = c(10:13), "DV-cells" = 14, "AP-cells" = c(16:15),
  "amp-SCs" = c(7,9), "mantle-cells" = c(5,6), "skin" = c(3,4))

meta$cell.type.ident <- factor(rep("", nrow(meta)),
  levels = names(cells), ordered = TRUE)

for (i in 1:length(cells)) {
  meta$cell.type.ident[meta$tree.ident %in% cells[[i]]] <- names(cells)[i]
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


p4 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.20,
  label = TRUE, label.size = 4, group.by = "cell.type.old")
p4 <- cleanUMAP(p4)

png(figurePath(paste(
  "anchored_cell_type_old.png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p4)
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
  saveRDS(obj_integrated, dataPath(paste0(
    "SeurObj_anchored_cell_type_update_", script_name,"_.RDS")))
  obj_integrated <- readRDS(dataPath(paste0(
    "SeurObj_anchored_cell_type_update_", script_name,"_.RDS")))
}

str(obj_integrated@meta.data)
dim(obj_integrated@assays$RNA@data)
grep("^sox2$", rownames(obj_integrated), value = TRUE)


# ========================================================= Additional figures
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
DefaultAssay(obj_integrated) <- "RNA"

Idents(obj_integrated) <- obj_integrated@meta.data$integrated_snn_res.1.8
levels(Idents(obj_integrated))
obj_integrated <- BuildClusterTree(
  obj_integrated, reorder = TRUE, reorder.numeric = TRUE )


# ==== mclapply for Find Markers
marker_table <- mcFindMarkers(obj_integrated,
  save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "cluster_dim10_res0.6_")

WriteXLS::WriteXLS(marker_table,
  figurePath(paste0("cluster_markers_dim10_res0.6_", script_name,".xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)


# ==== mclapply for plotting cluster markers
mcPlotMarkers(obj_integrated, marker_table, n_genes = 200,
  folder_prefix = "clusters_dim10_res0.6_")


# ==== treatments by diff exp
Idents(obj_integrated) <- obj_integrated@meta.data$cell.type.ident
levels(Idents(obj_integrated))

all_cells <- diffConditionMrkrs(obj_integrated, n_cores = 10)

WriteXLS::WriteXLS(all_cells,
  figurePath(paste0("all_cell_types_compared_all_cond.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

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
all_cells <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "all_cells_condition_diff_nmast_skin_additional_timepoints_v1.0_.RDS"))

all_cells <- all_cells[order(all_cells$p_val_adj,
  all_cells$cell.type.and.trt),]

diffConditionPlots(obj_integrated, input_file = all_cells,
  n_cores = 10, n_genes = 400, folder_prefix = "all_cells")

central_markers <- readRDS(
  paste0(""))

central_markers <- central_markers[order(central_markers$p_val_adj,
  central_markers$cell.type.and.trt),]

diffConditionPlots(obj_integrated, input_file = central_markers,
  n_cores = 10, n_genes = 400, folder_prefix = "central-all-timepoints")

central_prog <- readRDS(
  paste0(""))

central_prog <- central_prog[order(central_prog$p_val_adj,
  central_prog$cell.type.and.trt),]

diffConditionPlots(obj_integrated, central_prog, n_cores = 10, n_genes = 400)

pbulk <- readRDS(
  paste0(""))

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


# ==== Percent cell type by treatment
meta <- obj_integrated@meta.data
table(meta$tree.ident)

for (i in seq_along(ids)) {
  print(ids[i])
  print(table(meta$tree.ident == 22 & meta$data.set == ids[i])[2] / 
    table(meta$data.set)[ids[i]])
  cat("\n")
}


# ==== 3D plot
# Create 3D object first in R-env on server ()

dims <- 1:10
obj_3D <- RunUMAP(obj_integrated, n.components = 3, dims = dims,
  min.dist = 0.3, spread = 1)
Idents(obj_3D) <- obj_3D@meta.data$cell.type.ident

levels(Idents(obj_3D))
obj_3D <- BuildClusterTree(
  obj_3D, reorder = TRUE, reorder.numeric = TRUE)

saveRDS(obj_3D, dataPath(paste0("SeurObj_anchored_3D_", script_name,".RDS")))

# (run locally)
# Load locally for plotly

# obj_3D <- readRDS(paste0("/Volumes/projects/ddiaz/Analysis/Data/sb2191-regen/",
#   "SeurObj_anchored_3D_nmast_skin_additional_timepoints_seurat3_v1.0.RDS"))

library(plotly)
library(dplyr)
library(Seurat)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Idents(obj_3D) <- obj_3D@meta.data$integrated_snn_res.1.8
obj_3D <- BuildClusterTree(
  obj_3D, reorder = TRUE, reorder.numeric = TRUE)

levels(Idents(obj_3D))
ident <- "cell.type.ident"

umap_dims <- list()
for (i in 1:3) {umap_dims[[i]] <-  obj_3D[["umap"]]@cell.embeddings[,i]}

plot_data <- FetchData(object = obj_3D,
  vars = c("UMAP_1", "UMAP_2", "UMAP_3", ident))

cell_colors <- gg_color_hue(max(as.numeric(plot_data[,ident])))

if (FALSE) {
  plot_ly(data = plot_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    mode = "markers", marker = list(size = 2, width = 2),
    color = ~cell.type.ident, text=~ident, colors = cell_colors,
    type = "scatter3d") %>%  hide_colorbar()
}

p <- plot_ly(data = plot_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
  mode = "markers", marker = list(size = 2, width = 2),
  color = ~cell.type.ident, text=~cell.type.ident, colors = cell_colors,
  type = "scatter3d") %>%  hide_colorbar()

htmlwidgets::saveWidget(as_widget(p), figurePath("3D_UMAP_res1.80_.html"))