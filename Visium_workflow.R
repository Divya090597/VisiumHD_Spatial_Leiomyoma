options(repos = c(CRAN = "https://cloud.r-project.org/"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


# VisiumHD_workflow -------------------------------------------------------

seu_obj = readRDS("/home/rstudio/VisiumHD_LEO/LEO_Annotation_Updated_Seurat_2142025.rds")
VlnPlot(obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt"), ncol = 3)

summary(obj@meta.data[c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt")])

obj = subset(obj, subset = nFeature_Spatial.008um > 50 & obj$nCount_Spatial.008um > 0 & !is.na(obj@meta.data[["percent.mt"]]))

obj = subset(obj, subset = nFeature_Spatial.008um > 0 & nCount_Spatial.008um > 0 & !is.na(obj@meta.data[["percent.mt"]]))

obj@images  

VlnPlot(obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt"), ncol = 3)

obj = subset(obj, subset = percent.mt < 50)

SpatialFeaturePlot(seu_obj, features = c("S1OOA6","S100A9"))


# Plot the two genes separately and combine them in one layout
p1 <- SpatialFeaturePlot(
  object = seu_obj,
  features = "CD3D",                # First gene
  pt.size.factor = 5,               # Adjust point size
  alpha = c(0.8, 1),                # Transparency for low and high expression
  image.alpha = 0.5                 # Reduce background image opacity
) +
  scale_fill_gradient(low = "grey", high = "red") +
  labs(title = "Spatial Expression of CD3") +
  theme_minimal()

p2 <- SpatialFeaturePlot(
  object = seu_obj,
  features = "MS4A1",                # Second gene
  pt.size.factor = 5,               # Adjust point size
  alpha = c(0.8, 1),                # Transparency for low and high expression
  image.alpha = 0.5                 # Reduce background image opacity
) +
  scale_fill_gradient(low = "grey", high = "blue") +
  labs(title = "Spatial Expression of Bcell ") +
  theme_minimal()
p2
# Combine the two plots side-by-side
library(patchwork)
p1 + p2

SpatialDimPlot(
  object = seu_obj,
  group.by = "Cell_Type",  # Replace with your grouping column
  label = TRUE,                  # Show cluster labels
  label.size = 3,                # Adjust label size
  pt.size.factor = 3,          # Increase point size for clarity
  alpha = c(0.8, 1),             # Make points slightly transparent for better overlay
  image.alpha = 0.5,             # Reduce image opacity for better cluster visibility
) +
  labs(title = "Spatial DimPlot") +
  theme_minimal()

SpatialDimPlot(
  object = seu_obj,
  group.by = "Cell_Type",  # Replace with your grouping column
  label = TRUE,                  # Show cluster labels
  label.size = 3,                # Adjust label size
  pt.size.factor = 3,            # Increase point size for clarity
  alpha = c(0.8, 1),             # Make points slightly transparent for better overlay
  image.alpha = 0.5              # Reduce image opacity for better cluster visibility
) +
  labs(title = "Spatial DimPlot") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),  # Bold title
    axis.title = element_text(face = "bold", size = 12),  # Bold axis labels
    axis.text = element_text(face = "bold", size = 10),   # Bold axis text
    legend.text = element_text(face = "bold", size = 10), # Bold legend text
    legend.title = element_text(face = "bold", size = 12) # Bold legend title
  )

Idents(seu_obj) = "Cell_Type"
SpatialDimPlot(
  object = seu_obj,
  cells.highlight = WhichCells(seu_obj, idents = c("Tcells")),  
  cols.highlight = c("red", "grey"),
  pt.size.factor = 3,
  image.alpha = 0.5
)

SpatialDimPlot(
  object = seu_obj,
  cells.highlight = list(
    "Fibroblasts" = WhichCells(seu_obj, idents = "Fibroblasts"),
    "Endothelialcells" = WhichCells(seu_obj, idents = "Endothelialcells")
  ),
  cols.highlight = c("red", "blue","gray"),  # Add a color for non-highlighted cells
  pt.size.factor = 3,                        # Adjust point size for visibility
  image.alpha = 0.5                          # Reduce image opacity to emphasize the highlights
) +
  labs(title = "Spatial DimPlot: Fibroblasts vs Tcells") +
  theme_minimal()


SpatialDimPlot(
  object = seu_obj,
  cells.highlight = list(
    "SMC_a" = WhichCells(seu_obj, idents = "Smoothmuscle_a"),
    "SMC_b" = WhichCells(seu_obj, idents = "Smoothmuscle_b")
  ),
  cols.highlight = c("red", "blue","gray"),  # Add a color for non-highlighted cells
  pt.size.factor = 3,                        # Adjust point size for visibility
  image.alpha = 0.5                          # Reduce image opacity to emphasize the highlights
) +
  labs(title = "Spatial DimPlot: Smoothmuscle_a vs Smoothmuscle_b") +
  theme_minimal()

Idents(seu_obj) = "Cell_Type"
SpatialDimPlot(
  object = seu_obj,
  cells.highlight = WhichCells(seu_obj, idents = c("1_Trans_COLhigh", "36_Renewing","31_Renewing","41_Transitioning" ,"20_Renewing", "19_Renewing", "22_Transitioning","12_Renewing","11_SelfRenewing", "30_Transitioning","18_Renewing","17_Renewing", "2_Transitioning" , "7_Transitioning", "24_Transitioning","42_Transition","14_Renewing", "26_Renewing", "21_Renewing","23_Transitioning","39_Transition", "40_Renewal","25_Transitioning","31_Renewing","33_Transitioning","35_Transitioning","34_Transitioning","45_Transition","13_Renewing")),
  cols.highlight = c( "red","grey"),  # Highlight selected clusters in red
  pt.size.factor = 2,
  image.alpha = 0.5,
  label = T,
  label.size = 4
) +
  labs(title = "Spatial DimPlot Highlighting Tumor Clusters") +
  theme_minimal()

obj = NormalizeData(obj)
obj = FindVariableFeatures(obj)
obj = ScaleData(obj)
obj = RunPCA(obj)
options(future.globals.maxSize = 2000 * 1024^2)  # Set to 2 GB
obj = FindNeighbors(obj, dims = 1:50)
obj = FindClusters(obj, resolution = 0.5)
obj = RunUMAP(obj, dims = 1:30)
DimPlot(obj, reduction = "umap")

saveRDS(seu_obj, file = "/home/rstudio/VisiumHD_LEO/LEO_Annotation_Updated_Seurat.rds")

object@meta.data$sketch_snn_res.0.3_annotation <- NULL


object <- SketchData(
  object = obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(object) <- "sketch"

object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched_0.3", resolution = 0.3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
DimPlot(seu_obj, reduction = "umap.sketch", group.by = "Rough_Annotation", label = T)+
  theme(
    axis.title = element_text(face = "bold",size = 12),   
    axis.text = element_text(face = "bold",size = 12),    
    legend.text = element_text(face = "bold",size = 12),  
    legend.title = element_text(face = "bold",size = 12)  
  )

library(future)

# Increase global size limit
options(future.globals.maxSize = 2 * 1024^3)  # Increase to 2 GB

# Project the data
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched_0.3")
)


SpatialDimPlot(object, label = T)

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched_0.3"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(seu_obj) <- "Spatial.008um"
Idents(seu_obj) <- "seurat_cluster.projected"
p2 <- DimPlot(seu_obj, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2
Idents(seu_obj) <- "seurat_cluster.projected"
cells <- CellsByIdentities(seu_obj, idents = c(3))
p <- SpatialDimPlot(seu_obj,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
seu_obj2 <- BuildClusterTree(seu_obj, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p
Idents(object) = "seurat_clusters" 
library(Matrix)

agg <- AggregateExpression(object,
                                group.by =c( "orig.ident","Spatial.008um_snn_res.0.5",""),
                                assays = "sketch",
                                slot = "counts",
                                return.seurat = T)

View(agg@meta.data)
myofibro = c("MYH11","MYL9","TAGLN","TPM1","TPM2","ACTG2","ACTA2","TNS1","PALLD","CNN1","DES","CALD1","CAPN6","HMGA2","SMOC1","IGF2BP3","IGF2BP2","IL17B","MYH3","ESR1","PGR","VIM","LUM","ALDH1","CD90","FN1","DCN","OGN","MGP","COL1A1","COL1A2","COL3A1","PECAM1","CDH5","VWF","VCAM1")
Mono = c("ITGAM","ITGAX","CCR2","CD14","FCGR3A","VSIR","S100A8","S100A9","TLR2","PECAM1","CD163","CD86","MRC1","C1QA","C1QB","CD68","SPP1","HLA-DRA","HLA-DRB1")
Tumor = c("CAPN6","HMGA2","SMOC1","IGF2BP3","IGF2BP2","IL17B","MYH3","ESR1","PGR","CEA","MUC16","SCC")
gen = c("PTPRC","CD3D","CD3E","ADH1A","CD4","CD8A","CTLA4","FOXP3","ITGAM","ITGAX","CCR2","CD14","S100A8","S100A9","TLR2","CD163","FCGR3A","CD86","MRC1","C1QA","C1QB","CD68","SPP1","IRF8","IRF4","CD161","NKG7","GNLY","GZMB","GZMA","PRF1","CST7","CD7","DPP4","FCGR1A","HLA-DRA","HLA-DRB1")
stroma = c("MYC","MKI67","CD44","SOX2","C3","C7","S100A10","SCARA5","SFRP2","PI16","VIM","LUM","DCN","PDPN","VCAN","FN1","COL1A1","COL3A1","COL1A2","COL4A5","COL4A6","COL6A1","COL6A2","COL6A3","PECAM1","VWF","CDH5","VCAM1","TAGLN","ACTA2","CNN1","DES","TPM2","CAPN6","MED12","HMGA2","FH","SMOC1","IGF2BP3","IGF2BP2","IL17B","MYH3","TGFB1","TGFB3","ACVR1","ESR1","PGR","ACP1","PTH")
genes = c("PTPRC","EOMES","FOXP3","CD3D","CD3E","CD4","CD8A","CD19","IGHA1","MS4A1","ITGAM","ITGAX","CCR2","CD14","FCGR3A","S100A8","S100A9","TLR2","CD163","MRC1","SPP1","C1QB","IRF8","IRF4","XCR1","CLEC10A","NKG7","GNLY")

lipid = c("ALOX12","ALOX15B","ALOX5","ALOX5AP","LTA4H","PTGDS","PTGER1","PTGER2","PTGER3","PTGER4","PDGFA","EGF","CSF1")

DotPlot(seu_obj, features = lipid, group.by = "Cell_Type") +
  scale_size(range = c(3, 10)) + 
  scale_color_gradient(low = "white", high = "red") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  
    axis.text.y = element_text(size = 14, face = "bold")
  )

DefaultAssay(seu_obj) = "Spatial.008um"
# Extract scaled expression data for the selected genes
expression_data <- GetAssayData(seu_obj, slot = "scale.data")

col_fun <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))

cell_types <- seu_obj@meta.data$Annotation_2




Idents(seu_obj) = "seurat_cluster.projected"
seu_obj = RenameIdents(seu_obj, '0'= "Endothelialcells",'1' = "Artifact",'2'="Tumor" ,'3'="TAMS",'4'="Smoothmuscle_a", '5'="CD31high_Monocytes", '6'="Smoothmuscle_a",'7'="Tumor",'8'="EC_Programing",'9'="EC_Programing",'10'="Fibroblasts",'11'= "KI67high_Tumor",'12'="Tumor",'13'="Tumor",'14'="Tumor",'15'= "Smoothmuscle_a",'16'="Fibroblasts","17" = "Tumor","18"="Tumor","19"="Tumor","20"="Tumor","21"="Tumor","22"="Trans_Fibroblasts","23"="Tumor","24"="Tumor","25"="Trans_Fibroblasts","26"="Tumor","27"="Fibroblasts","28"="Tcells","29"="Smoothmuscle_a","30"="Trans_Fibroblasts","31"="Tumor","32"="Smoothmuscle_a","33"="Trans_Fibroblasts","34"="Tumor","35"="Tumor","36"="Tumor","37"="Trans_Fibroblasts","38"="Fibroblasts","39"="Trans_Fibroblasts","40"="Trans_Fibroblasts","41"="Tumor","42"="Tumor","43"="Endothelialcells","44"="Smoothmuscle_b","45"="Trans_Fibroblasts","46"="Endothelialcells","47"="Fibroblasts","48"="Fibroblasts" )
seu_obj@meta.data$Cell_Type = Idents(seu_obj)

Subset1 = subset(object, subset = seurat_cluster.sketched_0.3 %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29"))
Subset2 = subset(object, subset = seurat_cluster.sketched_0.3 %in% c("30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59"))

Sketched_markers_1 = FindMarkers(seu_obj, only.pos = T, min.pct = 0.001, logfc.threshold = 0.1)

fwrite(Sketched_markers,"/home/rstudio/VisiumHD_LEO/LEO_Markers_res.0.3.csv")


?FindAllMarkers

# Convert the factor column to character
seu_obj@meta.data$Annotation_2 <- as.character(seu_obj@meta.data$Annotation_2)
seu_obj@meta.data$Cell_Type = as.character(seu_obj@meta.data$Cell_Type)

# Update the value
seu_obj@meta.data$Annotation_2[seu_obj@meta.data$Annotation_2 == "1_Fibro_COLhigh"] <- "1_Artificial_pattern"
seu_obj@meta.data$Cell_Type[seu_obj@meta.data$Cell_Type == "Fibro_COLhigh"] = "Artificial_pattern"

# Convert back to factor if needed
seu_obj@meta.data$Annotation_2 <- as.factor(seu_obj@meta.data$Annotation_2)
seu_obj@meta.data$Cell_Type = as.factor(seu_obj@meta.data$Cell_Type)


# Replace NA values in the "Annotation_2" column with "1_Fibro_COLhigh"
seu_obj@meta.data$Annotation_2[is.na(seu_obj@meta.data$Annotation_2)] <- "1_Fibroblast_COLhigh"

# Subset the Seurat object to include only the specified cell types
selected_celltypes2 <- c("Tumor", "KI67high_Tumor")

Idents(seu_obj) = "Cell_Type"
# Create a subset of the Seurat object
subset <- subset(seu_obj, idents = selected_celltypes2)

# Spatial DimPlot for the subsetted object
SpatialDimPlot(
  object = subset,
  group.by = "Annotation_2",  # Adjust group.by to your column
  pt.size.factor = 4,                     # Adjust point size for better visibility
  image.alpha = 0.5,                      # Reduce histological image opacity
  label = TRUE,                           # Display labels for all clusters
  label.size = 4                          # Adjust label size for clarity
) +
  labs(title = "Spatial DimPlot of Tumor") +
  theme_minimal()

SpatialDimPlot(
  object = subset2,
  group.by = "Rough_Annotation",  # Use the metadata column for cluster annotations
  pt.size.factor = 3,             # Increase point size for better visibility
  image.alpha = 0.3,              # Reduce histological image opacity
  label = T,                   # Display labels for all clusters
  label.size = 5,                 # Make labels larger for clarity
  label.box = TRUE                # Add label boxes for better separation
) +
  labs(title = "Spatial DimPlot FOR Tumor") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),  # Bold and larger plot title
    legend.text = element_text(face = "bold", size = 12), # Bold legend text
    axis.text = element_text(face = "bold"),              # Bold axis text
    axis.title = element_text(face = "bold")              # Bold axis titles
  )

# Define myeloid clusters
myeloid_clusters <- c("5_CD31high_Monocytes", "3_Macrophage")

# Subset the Seurat object
myeloid_subset2 <- subset(seu_obj, subset = Rough_Annotation %in% myeloid_clusters)

# Normalize the subsetted data
myeloid_subset <- NormalizeData(myeloid_subset)

# Find variable features
myeloid_subset <- FindVariableFeatures(myeloid_subset)

# Scale the data
myeloid_subset <- ScaleData(myeloid_subset)

# Perform PCA
myeloid_subset <- RunPCA(myeloid_subset)

# Find neighbors and clusters (adjust resolution for granularity)
myeloid_subset <- FindNeighbors(myeloid_subset, dims = 1:30)
myeloid_subset <- FindClusters(myeloid_subset, resolution = 0.5)  # Adjust resolution as needed

# Run UMAP for visualization
myeloid_subset <- RunUMAP(myeloid_subset, dims = 1:30)

# Visualize the clusters
DimPlot(myeloid_subset, reduction = "umap",group.by = "Rough_Annotation" ,label = TRUE, label.size = 5) +
  labs(title = "UMAP of Reclustered Myeloid Subset")


