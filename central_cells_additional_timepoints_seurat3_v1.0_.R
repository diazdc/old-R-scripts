# 1.0 Subsetted cells from all_LL_cells_regen_anchored_seurat3_v1.2_.R

library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("central_cells_additional_timepoints_seurat3_v1.0")

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
    "SeurObj_anchored_cell_type_update_nmast_skin_",
    "additional_timepoints_seurat3_v1.0_.RDS"))

to_subset <- seurat_obj@meta.data$cell.type.ident %in% "central-cells"

seurat_obj <- seurat_obj[,to_subset]
colnames(seurat_obj@meta.data)[10] <- "cell.type.old"


# ================================================ Normalization/Dim Reduction


DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- NormalizeData(seurat_obj,
  normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- FindVariableFeatures(
  seurat_obj, selection.method = "vst", nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

pdf(figurePath("PC_elbowPlot.pdf"))
ElbowPlot(seurat_obj)
dev.off()

dims <- 1:10
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
  "anchored_clusters_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p1)
dev.off()


# ==== Plot by treatment
p2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.40,
  label = FALSE, label.size = 4, group.by = "data.set", cols = trt_colors)
p2 <- cleanUMAP(p2)

png(figurePath(paste("anchored_dataset_dims", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p2)
dev.off()

png(figurePath(paste0("UMI_count_by_dataset.png")),
  width = 12, height = 10, units = "in", res = 300)
VlnPlot(seurat_obj,"nCount_RNA",
  group.by= "data.set", pt.size = 0.3)
dev.off()

stop("End Analysis")

if(FALSE) {
  source(paste0("/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/",
    "Seurat3-final/nmast_skin_additional_timepoints_seurat3_v1.0_.R"),
    echo = TRUE)
}

if (FALSE){
  testPCs(seurat_obj, from_to = 5:9)
}


# ============================================================= Save/load data


if(FALSE) {
  saveRDS(seurat_obj,
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
  print("Object saved")

  seurat_obj <- readRDS(
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
}

table(seurat_obj@meta.data[
  seurat_obj@meta.data$data.set == "homeo", "tree.ident"])
table(seurat_obj@meta.data$tree.ident)

# check res first. using res 0.6 for cell classification
levels(Idents(seurat_obj))
colnames(seurat_obj@meta.data)
Idents(seurat_obj) <- seurat_obj@meta.data$integrated_snn_res.0.6

# Adding in cell type names
meta <- seurat_obj@meta.data
colnames(meta)

cells <- list("mature-HCs" = 1, "early-HCs" = 2,  "HC-prog" = 6,
  "central-cells" = c(10:13), "DV-cells" = 14, "AP-cells" = c(16:15),
  "amp-SCs" = c(5,7), "mantle-cells" = c(3,4), "col1a1b-pos" = c(9),
  "c1qtnf5-pos" = 8, "clec14a-pos" = 10)

meta$cell.type.ident <- factor(rep("", nrow(meta)),
  levels = names(cells), ordered = TRUE)

for (i in 1:length(cells)) {
  meta$cell.type.ident[meta$tree.ident %in% cells[[i]]] <- names(cells)[i]
}

seurat_obj@meta.data <- meta
Idents(seurat_obj) <- seurat_obj@meta.data$cell.type.ident

p3 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.20,
  label = TRUE, label.size = 4)
p3 <- cleanUMAP(p3)

png(figurePath(paste(
  "anchored_cell_type_res", res, "dim", max(dims),".png", sep = "_")),
  width = 12, height = 10, units = "in", res = 300)
print(p3)
dev.off()


# ==== median gene count per cluster
vln_list <- list()[1:9]
vln_list <- lapply(seq_along(ids), function(i) {
  ind <- which(seurat_obj@meta.data$data.set == ids[i])
  VlnPlot(seurat_obj[,ind],"nCount_RNA",
    group.by= "tree.ident", pt.size = 0.3) + ggtitle(ids[i])
  })

png(figurePath(paste0("UMI_count_by_cluster.png")),
  width = 25, height = 15, units = "in", res = 300)
print(cowplot::plot_grid(plotlist = vln_list), ncol = 3)
dev.off()

# ! Final Object
if(TRUE) {
  saveRDS(seurat_obj, dataPath(paste0(
    "SeurObj_anchored_cell_type_update", script_name,"_.RDS")))
  seurat_obj <- readRDS(dataPath(paste0(
    "SeurObj_anchored_cell_type_update", script_name,"_.RDS")))
}


# ========================================================= Additional figures
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
DefaultAssay(seurat_obj) <- "RNA"

Idents(seurat_obj) <- seurat_obj@meta.data$integrated_snn_res.1.8
levels(Idents(seurat_obj))
seurat_obj <- BuildClusterTree(
  seurat_obj, reorder = TRUE, reorder.numeric = TRUE )


# ==== mclapply for Find Markers
marker_table <- mcFindMarkers(seurat_obj,
  save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "cluster_dim10_res0.6_")

WriteXLS::WriteXLS(marker_table,
  figurePath(paste0("cluster_markers_dim10_res0.6_", script_name,".xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)


# ==== mclapply for plotting cluster markers
mcPlotMarkers(seurat_obj, marker_table, n_genes = 200,
  folder_prefix = "clusters_dim10_res0.6_")


# ==== treatments by diff exp
Idents(seurat_obj) <- seurat_obj@meta.data$cell.type.ident
levels(Idents(seurat_obj))

all_cells <- diffConditionMrkrs(seurat_obj, n_cores = 10)

WriteXLS::WriteXLS(all_cells,
  figurePath(paste0("all_cell_types_compared_all_cond.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

single_group <- diffConditionMrkrs(seurat_obj, n_cores = 1,
  group_clusters = c("central-cells"), file_prefix = "central_cells")

WriteXLS::WriteXLS(single_group,
  figurePath(paste0("AP_combined_vs_other_all_timepoints.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

mantle_amp_AP <- diffConditionMrkrs(seurat_obj, n_cores = 1,
  group_clusters = c("mantle-cells", "amp-SCs", "AP-cells"))

WriteXLS::WriteXLS(mantle_amp_AP,
  figurePath(paste0("mantle_amp_AP_combined_vs_other_all_timepoints.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

pbulk <- diffConditionMrkrs(seurat_obj, n_cores = 1,
  cell_group_name = "pbulk", group_clusters = levels(Idents(seurat_obj)))

WriteXLS::WriteXLS(pbulk,
  figurePath(paste0("pbulk_all_timepoints.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)


# ==== plot treatment by diff exp
all_cells <- readRDS(
  paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "all_cells_condition_diff_nmast_skin_additional_timepoints_v1.0_.RDS"))

all_cells <- all_cells[order(all_cells$p_val_adj,
  all_cells$cell.type.and.trt),]

diffConditionPlots(seurat_obj, input_file = all_cells,
  n_cores = 10, n_genes = 400, folder_prefix = "all_cells")

central_markers <- readRDS(
  paste0(""))

central_markers <- central_markers[order(central_markers$p_val_adj,
  central_markers$cell.type.and.trt),]

diffConditionPlots(seurat_obj, input_file = central_markers,
  n_cores = 10, n_genes = 400, folder_prefix = "central-all-timepoints")

central_prog <- readRDS(
  paste0(""))

central_prog <- central_prog[order(central_prog$p_val_adj,
  central_prog$cell.type.and.trt),]

diffConditionPlots(seurat_obj, central_prog, n_cores = 10, n_genes = 400)

pbulk <- readRDS(
  paste0(""))

pbulk <- pbulk[order(pbulk$p_val_adj, pbulk$cell.type.and.trt),]

diffConditionPlots(seurat_obj, pbulk, all_idents = TRUE,
  folder_prefix = "pbulk", n_cores = 10, n_genes = 400)


# ==== Correlated genes
cor_list <- list()

sox4a_cor_list <- lapply(1:9, function(i) {
    cor_list[[i]] <- corrleatedGenes(seurat_obj, "sox4a*1",
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
  DefaultAssay(seurat_obj) <- "RNA"
  assay <- "assay_RNA_"
}

seurat_obj <- ScaleData(seurat_obj)
all_markers <- readRDS(dataPath(paste0("AllMarkers_", script_name,"_.RDS")))

# Subset by identity
meta <- seurat_obj@meta.data

treatment <- "central-cells"
BCs_to_subset <- rownames(meta[meta$cell.type.ident == treatment,])

# Subset object
obj_subset <- seurat_obj[,BCs_to_subset]

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
DefaultAssay(seurat_obj) <- "RNA"

common_features <- scan(paste0("/n/projects/ddiaz/Analysis/",
  "Data/gene-lists/common_features.txt"), what = "character")

e <- FeaturePlot(seurat_obj, common_features,
  reduction = "umap", pt.size = 0.25, combine = FALSE)

for (i in 1:length(e)) {
  e[[i]] <- e[[i]] + NoLegend() + NoAxes()
}

png(figurePath("common_features.png"), width = 40,
  height = 80, units = "in", res = 200)
print(cowplot::plot_grid(plotlist = e, ncol = 4))
dev.off()


# ==== Percent cell type by treatment
meta <- seurat_obj@meta.data
table(meta$tree.ident)

for (i in seq_along(ids)) {
  print(ids[i])
  print(table(meta$tree.ident == 22 & meta$data.set == ids[i])[2] / 
    table(meta$data.set)[ids[i]])
  cat("\n")
}


# ==== 3D plot
# Create 3D object first in R-env on server ()

dims <- 1:11
obj_3D <- RunUMAP(seurat_obj, n.components = 3, dims = dims,
  min.dist = 0.3, spread = 1)
Idents(obj_3D) <- obj_3D@meta.data$integrated_snn_res.1.8

levels(Idents(obj_3D))
obj_3D <- BuildClusterTree(
  obj_3D, reorder = TRUE, reorder.numeric = TRUE)

saveRDS(obj_3D, dataPath(paste0("SeurObj_anchored_3D_", script_name,".RDS")))

# (run locally)
# Load locally for plotly
obj_3D <- readRDS(paste0("/Volumes/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "SeurObj_anchored_3D_nmast_skin_additional_timepoints_seurat3_v1.0.RDS"))

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
ident <- "tree.ident"

umap_dims <- list()
for (i in 1:3) {umap_dims[[i]] <-  obj_3D[["umap"]]@cell.embeddings[,i]}

plot_data <- FetchData(object = obj_3D,
  vars = c("UMAP_1", "UMAP_2", "UMAP_3", ident))

cell_colors <- gg_color_hue(max(as.numeric(plot_data[,ident])))

plot_ly(data = plot_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
  mode = "markers", marker = list(size = 2, width = 2),
  color = ~tree.ident, text=~ident, colors = cell_colors,
  type = "scatter3d") %>%  hide_colorbar()

p <- plot_ly(data = plot_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
  mode = "markers", marker = list(size = 2, width = 2),
  color = ~tree.ident, text=~ident, colors = cell_colors,
  type = "scatter3d") %>%  hide_colorbar()

htmlwidgets::saveWidget(as_widget(p), figurePath("3D_UMAP_res1.80_.html"))

sample_genes <- c("ahnak", "prkd2", "si:cabz01007794.1", "ptp4a2b", "ef1", "rhoca", "FO904881.1")

seurat_obj <- ScaleData(seurat_obj)
table(sample_genes %in% rownames(seurat_obj))

DoHeatmap(seurat_obj, "sample_genes")


# ==== biomart for human gene names
library(biomaRt)

test <- marker_table$Gene.stable.ID

mart = useMart(host = "http://sep2019.archive.ensembl.org", 
  biomart = "ensembl", dataset = "drerio_gene_ensembl")

mart_results <- getBM( values = test, filters = 'ensembl_gene_id',
  attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name"),
  mart = mart, uniqueRows = FALSE)

listAttributes(mart)
