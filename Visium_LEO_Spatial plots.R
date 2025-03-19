library(Seurat)
object=readRDS("~/VisiumHD/Seu_obj_updated.rds")

DimPlot(obj, reduction = "umap",label = T, pt.size = 0.5)
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters")
head(obj@reductions$umap@cell.embeddings)
obj <- RunUMAP(obj, dims = 1:30)
View(obj@meta.data)

obj@reductions$umap <- NULL


dir.create("C:/Users/80029518/Documents/VisiumHD/Leo2_625M/Leo2_625M/spatial", recursive = TRUE)
untar("C:/Users/80029518/Documents/VisiumHD/Leo2_625M/Leo2_625M/spatial.tar.gz",
      exdir = "C:/Users/80029518/Documents/VisiumHD/Leo2_625M/Leo2_625M/spatial")


# Extract
untar("~/VisiumHD/Leo2_625M/Leo2_625M/binned_outputs.tar.gz",
      exdir = "~/VisiumHD/Leo2_625M/Leo2_625M/binned_outputs")



# Load libraries ----------------------------------------------------------
library(Seurat)
library(Matrix)
library(hdf5r)
library(rhdf5)
library(arrow)
library(ggplot2)
library(patchwork)
library(dplyr)

data.dir = "~/VisiumHD/Leo2_625M/Leo2_625M"
object = Load10X_Spatial(data.dir = data.dir, bin.size = 8)

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
View(object@meta.data)

SpatialFeaturePlot(object, features = "nCount_Spatial.008um")
VlnPlot(object, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt"), ncol = 3)
object <- subset(object, subset = nFeature_Spatial.008um > 500 & nFeature_Spatial.008um < 4000 & nCount_Spatial.008um < 20000 & percent.mt < 10)


saveRDS(object,"~/VisiumHD/Visum_Seurat_Object.rds")



obj = readRDS("~/VisiumHD/Visum_Seurat_Object.rds")
VlnPlot(obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt"), ncol = 3)

summary(obj@meta.data[c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt")])

obj = subset(obj, subset = nFeature_Spatial.008um > 10 & obj$nCount_Spatial.008um > 0 & !is.na(obj@meta.data[["percent.mt"]]))

obj = subset(obj, subset = nFeature_Spatial.008um > 0 & nCount_Spatial.008um > 0 & !is.na(obj@meta.data[["percent.mt"]]))

obj@images  

VlnPlot(obj, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt"), ncol = 3)

obj = subset(obj, subset = percent.mt < 20)

SpatialFeaturePlot(obj, features = "nCount_Spatial.008um")
p1 = SpatialFeaturePlot(object, features = c("CD14","CD163","S100A3"))
p1
p2 = SpatialFeaturePlot(object, features = c("TNC","PECAM1"))
p2
p3 = SpatialFeaturePlot(object, features = c("FN1","VEGFA","COL3A1"))
p3
p4 = SpatialFeaturePlot(object, features = c("FABP4","PPARG"))
p4
SpatialDimPlot(object, label = T, repel = T, label.size = 4)
