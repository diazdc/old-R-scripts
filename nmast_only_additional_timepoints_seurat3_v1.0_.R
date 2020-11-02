sw# 1.0 Subsetted cells from anchored_cell_type_updateadditional_timepoints_seurat3_v1.0

{
library(Seurat)
library(WriteXLS)
library(ggplot2)
library(dplyr)
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

# Clear workspace
rm(list = ls())

script_name <- paste0("nmast_only_additional_timepoints_seurat3_v1.0")

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
  levels(Idents(obj_integrated))[c(1:8)]

obj_integrated <- obj_integrated[,to_subset]
colnames(obj_integrated@meta.data)[10] <- "cell.type.old"

table(obj_integrated@meta.data$nCount_RNA < 12000)
obj_integrated <- obj_integrated[,obj_integrated@meta.data$nCount_RNA < 12000]


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

if (FALSE) {
  DefaultAssay(obj_integrated) <- "integrated"
}

dims <- 1:10
res <- 0.9
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
    "Seurat3-final/nmast_only_additional_timepoints_seurat3_v1.0_.R"),
    echo = TRUE)
}

if (FALSE) {
  testPCs(obj_integrated, from_to = 16:25)
}


# ============================================================= Save/load data


if(FALSE) {
  saveRDS(obj_integrated,
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
  print("Object saved")

  obj_integrated <- readRDS(
    dataPath(paste0("SeurObj_", script_name,"_.RDS")))
}

# check res first. using res 0.9 for cell classification
levels(Idents(obj_integrated))
colnames(obj_integrated@meta.data)

# Adding in cell type names
meta <- obj_integrated@meta.data
colnames(meta)

cells <- list("mature-HCs" = 1, "early-HCs" = c(2,3),  "HC-prog" = c(4,9),
  "central-cells" = c(11:15,19), "DV-cells" = 17, "AP-cells" = c(16,20),
  "amp-SCs" = c(8,10,18), "mantle-cells" = c(5:7))

meta$cell.type.ident <- factor(rep("", nrow(meta)),
  levels = names(cells), ordered = TRUE)

for (i in 1:length(cells)) {
  meta$cell.type.ident[meta$tree.ident %in% cells[[i]]] <- names(cells)[i]
}

obj_integrated@meta.data <- meta
Idents(obj_integrated) <- obj_integrated@meta.data$cell.type.ident

p3 <- DimPlot(obj_integrated, reduction = "umap", pt.size = 0.40,
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

str(obj_integrated)
dim(obj_integrated@assays$RNA@data)
grep("sox", rownames(obj_integrated), value = TRUE)


# Subset cell types
if (FALSE) {
  ap_obj <- obj_integrated[,
    obj_integrated@meta.data$cell.type.ident == "AP-cells"]

  dv_obj <- obj_integrated[,
    obj_integrated@meta.data$cell.type.ident == "DV-cells"]

  mantle_obj <- obj_integrated[,
      obj_integrated@meta.data$cell.type.ident == "mantle-cells"]
  
  sub_type_list <- list(AP = ap_obj, DV = dv_obj, mantle = mantle_obj)

  saveRDS(sub_type_list, dataPath(paste0(
    "SeurObj_list_AP_DV_mantle", script_name,"_.RDS")))

  dv_obj <- obj_integrated[,
    obj_integrated@meta.data$cell.type.ident == "DV-cells"]
  saveRDS(dv_obj, dataPath(paste0(
    "SeurObj_DV_cells_", script_name,"_.RDS")))

  mantle_obj <- obj_integrated[,
      obj_integrated@meta.data$cell.type.ident == "mantle-cells"]
  saveRDS(mantle_obj, dataPath(paste0(
    "SeurObj_mantle_cells_", script_name,"_.RDS")))
}


# ========================================================= Additional figures
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
DefaultAssay(obj_integrated) <- "RNA"

Idents(obj_integrated) <- obj_integrated@meta.data$integrated_snn_res.1.8
levels(Idents(obj_integrated))


# ==== mclapply for Find Markers
marker_table <- mcFindMarkers(obj_integrated,
  save_raw = TRUE, pval_cutoff = 0.05, file_prefix = "cluster_dim10_res0.9_")

WriteXLS::WriteXLS(marker_table,
  figurePath(paste0("cluster_markers_dim10_res0.9_", script_name,".xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

mcPlotMarkers(obj_integrated, marker_table, n_genes = 200,
  folder_prefix = "clusters_dim10_res0.9_")


# ==== Compare conditions
devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")

DefaultAssay(obj_integrated) <- "RNA"
all_markers <- diffConditionMrkrs(obj_integrated, n_cores = 10)

WriteXLS::WriteXLS(all_markers,
  figurePath(paste0("timepoint_cell_type_markers.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

all_markers <- readRDS(paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "all_markers_nmast_only_additional_timepoints_seurat3_v1.0_.RDS"))
dim(all_markers)

diffConditionPlots(obj_integrated, input_table = all_markers,
  n_genes = 200, n_cores = 10, folder_prefix = "")

{# Creating multiple spreadsheets with separate tabs for easier veiwing
library(openxlsx)
split_df <- split(all_markers, all_markers$cell.type.and.trt)
cell_types <- unique(obj_integrated@meta.data$cell.type.ident)
colnames(all_markers)

diff_cond_dir <- figurePath("timepoint-cell-type-markers/")
dir.create(diff_cond_dir)

parallel::mclapply(seq_along(cell_types), function(h){
  wb <- createWorkbook()
  split_df_sub <- split_df[grep(cell_types[h], names(split_df))]
  
  # Trick for selecting every other element
  to_order <- unlist(strsplit(names(split_df_sub), "_"))[c(FALSE, TRUE)]
  split_df_sub <- split_df_sub[match(ids, to_order)]

  for (i in seq_along(split_df_sub)) {
    colnames(split_df_sub[[i]]) <- colnames(all_markers)
    print(colnames(split_df_sub[[i]]))
    addWorksheet(wb, names(split_df_sub[i]))
    writeData(wb, names(split_df_sub[i]), split_df_sub[[i]])
  }
  saveWorkbook(wb,
    paste0(diff_cond_dir, cell_types[h], ".xlsx"), overwrite = TRUE)
  }, mc.cores = 8
)
}

# Cell specific 
all_markers2 <- diffConditionMrkrs(
  obj_integrated, cell_specific = TRUE, n_cores = 10)

WriteXLS::WriteXLS(all_markers2,
  figurePath(paste0("timepoint_cell_specific_markers.xlsx")),
  row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

all_markers2 <- readRDS(paste0("/n/projects/ddiaz/Analysis/Data/sb2191-regen/",
  "all_markers_specific_nmast_only_additional_timepoints_seurat3_v1.0_.RDS"))


{# Creating multiple spreadsheets with separate tabs for easier veiwing
library(openxlsx)
split_df <- split(all_markers2, all_markers2$cell.type.and.trt)
cell_types <- unique(obj_integrated@meta.data$cell.type.ident)
colnames(all_markers2)

diff_cond_dir <- figurePath("timepoint-specific-markers/")
dir.create(diff_cond_dir)

parallel::mclapply(seq_along(cell_types), function(h){
  wb <- createWorkbook()
  split_df_sub <- split_df[grep(cell_types[h], names(split_df))]
  
  # Trick for selecting every other element
  to_order <- unlist(strsplit(names(split_df_sub), "_"))[c(FALSE, TRUE)]
  split_df_sub <- split_df_sub[match(ids, to_order)]

  for (i in seq_along(split_df_sub)) {
    colnames(split_df_sub[[i]]) <- colnames(all_markers2)
    print(colnames(split_df_sub[[i]]))
    addWorksheet(wb, names(split_df_sub[i]))
    writeData(wb, names(split_df_sub[i]), split_df_sub[[i]])
  }
  saveWorkbook(wb,
    paste0(diff_cond_dir, cell_types[h], ".xlsx"), overwrite = TRUE)
  }, mc.cores = 8
)
}

diffConditionPlots(obj_integrated, input_table = all_markers2,
  n_genes = 200, n_cores = 10, folder_prefix = "specific-")

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
# Load locally to launch plotly

# obj_3D <- readRDS(paste0("/Volumes/projects/ddiaz/Analysis/Data/sb2191-regen/",
#   "SeurObj_anchored_3D_nmast_skin_additional_timepoints_seurat3_v1.0.RDS"))
library(plotly)
library(dplyr)
library(Seurat)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

if (FALSE) {
  Idents(obj_3D) <- obj_3D@meta.data$integrated_snn_res.1.8
  obj_3D <- BuildClusterTree(
    obj_3D, reorder = TRUE, reorder.numeric = TRUE)
}

levels(Idents(obj_3D))
ident <- "cell.type.ident"

umap_dims <- list()
for (i in 1:3) {umap_dims[[i]] <- obj_3D[["umap"]]@cell.embeddings[,i]}

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

htmlwidgets::saveWidget(as_widget(p), figurePath("3D_UMAP_cell_type_.html"))

cell_sub <- rownames(obj_integrated@meta.data[
  obj_integrated@meta.data$data.set %in% "1hr" &
  obj_integrated@meta.data$cell.type.ident %in% "central-cells",])

FetchData(obj_integrated, vars = c("data.set"), cells = cell_sub)
VlnPlot(obj_integrated, features = "atoh1a", multi.group = cell_sub)



# ==== Stacked violing plot NOT by dataset
to_plot <- c("atoh1a", "atoh1b", "stmn1a", "pcna", "cldnb", "her4.1", "sost")

plot_list <- lapply(seq_along(to_plot), function(i) {
  p <- VlnPlot(obj_integrated, features = to_plot[i], pt.size = 0.0)  + 
  xlab("") + ylab(to_plot[i]) + ggtitle("") + 
  theme(legend.position = "none", axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = rel(1), angle = 0),
    axis.text.y = element_text(size = rel(1)),
    plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"))
  return(p)
  })

extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
  theme(axis.text.x=element_text(), axis.ticks.x = element_line())

ymaxs <- purrr::map_dbl(plot_list, extract_max)
ymaxs <- rep(max(ymaxs), length(ymaxs))
plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x +
  scale_y_continuous(breaks = c(y)) + expand_limits(y = y))

p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
p <- p + RotatedAxis()
​
png(figurePath("stackedviolin.png"), res = 200,
  units = "in", width = 10, height = 10)
print(p)
dev.off()


# ==== Stacked violin plot by DATASET
ids <- c("homeo", "15min", "1hr", "3hr", "5hr")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

obj_trt_list <- list()[1:5]
for (i in 1:5) {
  obj_trt_list[[i]] <- obj_integrated[,obj_integrated[["data.set"]] == ids[i]]
}

goi <- "pcna"
trt_plot_list <- list()[1:5]
names(trt_plot_list) <- ids
for (i in 1:5) {
  vln_obj <- VlnPlot(
    obj_trt_list[[i]], features = goi, pt.size = 0) +
    xlab("") + ylab(ids[i]) + ggtitle("") +
      theme(legend.position = "none", axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = rel(1), angle = 0),
      axis.text.y = element_text(size = rel(1)),
      plot.margin = unit(c(-0.75, 0.5, -0.75, 0.5), "cm"))
  trt_plot_list[[i]] <- vln_obj
}

extract_max <- function(p) {
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

trt_plot_list[[length(trt_plot_list)]] <-
  trt_plot_list[[length(trt_plot_list)]] +
  theme(axis.text.x = element_text(), axis.ticks.x = element_line())

# change the y-axis tick to only max value
ymaxs <- purrr::map_dbl(trt_plot_list, extract_max(vln_obj))
trt_plot_list <- purrr::map2(trt_plot_list, ymaxs, function(x, y) x +
  scale_y_continuous(breaks = c(y)) + expand_limits(y = y))

grid_obj <- cowplot::plot_grid(plotlist = trt_plot_list,
  nrow = 5, ncol = 1, axis = "l", align = "hv") +
  theme(plot.margin = margin(1.0, 1.0, 1.0, 1.0, unit = "in")) +
  ggtitle(goi) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

grid_obj2 <- cowplot::plot_grid(plotlist = trt_plot_list,
  nrow = 5, ncol = 1, axis = "l", align = "hv") +
  theme(plot.margin = margin(1.0, 1.0, 1.0, 1.0, unit = "in")) +
  ggtitle(goi) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# option for plotting more than
grid_obj3 <- cowplot::plot_grid(grid_obj, grid_obj2, ncol = 1)

png(figurePath(paste0(goi,"_stacked_vln.png")),
  width = 16, height = 8, units = "in", res = 300)
print(grid_obj2)
dev.off()

png(figurePath(paste0("test_stacked_vln.png")),
  width = 16, height = 16, units = "in", res = 300)
print(grid_obj3)
dev.off()


# ==== percentage of cells per timepoint
dataset_cnt <- table(obj_integrated@meta.data$data.set)

trts <- as.character(unique(obj_integrated@meta.data$data.set))
types <- sort(unique(obj_integrated@meta.data$cell.type.ident))

percents_by_treat <- lapply(seq_along(trts),
  function(h) {
    vapply(seq_along(types), function (i) {
      percentage <- table(obj_integrated@meta.data$data.set == trts[h] &
        obj_integrated@meta.data$cell.type.ident == types[i])[2] /
        dataset_cnt[trts[h]]
    return(percentage)}, double(1))
  }
)
names(percents_by_treat) <- ids
cell_types <- rep(types, 9)

percent_df <- plyr::ldply(percents_by_treat, data.frame)
colnames(percent_df) <- c("time", "percent")
percent_df <- cbind(cell_types, percent_df)
percent_df$time <- factor(percent_df$time, ordered = TRUE, levels  = ids)

# Stacked + percent
bar_plot <- ggplot(percent_df, aes(fill = cell_types, y = percent, x = time)) + 
    geom_bar(position = "fill", stat = "identity")

bar_plot <- bar_plot + theme_classic() +
  scale_fill_manual(values = gg_color_hue(8))
bar_plot <- bar_plot + guides(shape = guide_legend(order = 2))

pdf(figurePath("cell_proportion_barplot.pdf"), width = 8, height = 8)
print(bar_plot)
dev.off()

png(figurePath("cell_proportion_barplot.png"), width = 8, height = 8,
  units = "in", res = 300)
print(bar_plot)
dev.off()




