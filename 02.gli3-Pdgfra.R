##########Pdgfra###############

load("tdT_blank-scdata.Rdata")
load("tdT_sox-scdata.Rdata")
load("tdT_soxgli-scdata.Rdata")


setwd("E:\\single-cell\\data")
load("scdata1.Rdata")

library(clustree)
library(dplyr)
library(hdf5r)
library(Seurat)
library(data.table)
library(patchwork)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(harmony)

Convert("sorna_filtered.h5ad", dest = "h5Seurat", overwrite = TRUE)
seurat_object <- LoadH5Seurat("sorna_filtered.h5Seurat")



########
sgrna <- soxgli
sorna <- sox
blank <- blank

rm(tdT_positive_sgrna,tdT_positive_sorna,sgrna_data,sorna_data)

head(sgrna@meta.data)
head(sorna@meta.data)
sgrna[["nFeature_RNA"]] <- Matrix::colSums(sgrna@assays$RNA@counts > 0)
sgrna[["nCount_RNA"]] <- Matrix::colSums(sgrna@assays$RNA@counts)
sorna[["nFeature_RNA"]] <- Matrix::colSums(sorna@assays$RNA@counts > 0)
sorna[["nCount_RNA"]] <- Matrix::colSums(sorna@assays$RNA@counts)

blank[["nFeature_RNA"]] <- Matrix::colSums(blank@assays$RNA@counts > 0)
blank[["nCount_RNA"]] <- Matrix::colSums(blank@assays$RNA@counts)

metadata <- sorna@meta.data

sgrna[["orig.ident"]] <- "sgrna"
sorna[["orig.ident"]] <- "sorna"
blank[["orig.ident"]] <- "blank"

sgrna[["nFeature_RNA"]] <- Matrix::colSums(sgrna@assays$RNA@counts > 0)
sgrna[["nCount_RNA"]] <- Matrix::colSums(sgrna@assays$RNA@counts)
sorna[["nFeature_RNA"]] <- Matrix::colSums(sorna@assays$RNA@counts > 0)
sorna[["nCount_RNA"]] <- Matrix::colSums(sorna@assays$RNA@counts)
blank[["nFeature_RNA"]] <- Matrix::colSums(blank@assays$RNA@counts > 0)
blank[["nCount_RNA"]] <- Matrix::colSums(blank@assays$RNA@counts)

scRNA_pre_filter <- merge(x = sgrna, y = list(sorna, blank), 
                          add.cell.ids = c("sgrna", "sorna", "blank"))

scRNA_pre_filter[["nFeature_RNA"]] <- Matrix::colSums(scRNA_pre_filter@assays$RNA@counts > 0)
scRNA_pre_filter[["nCount_RNA"]] <- Matrix::colSums(scRNA_pre_filter@assays$RNA@counts)


scRNA_pre_filter[["mt_percent"]] <- PercentageFeatureSet(scRNA_pre_filter, pattern = "^mt-")
HB_genes_mouse <- c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-bh2", "Hbb-bh1")
HB_genes_mouse <- intersect(HB_genes_mouse, rownames(scRNA_pre_filter))
scRNA_pre_filter[["HB_percent"]] <- PercentageFeatureSet(scRNA_pre_filter, features = HB_genes_mouse)

VlnPlot(scRNA_pre_filter, 
              features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"), 
              pt.size = 0.01, 
              ncol = 4, 
              group.by = "orig.ident")

theme(axis.title.x = element_text(face="bold")) + 
  scale_x_discrete(labels = c("sgrna" = "sgrna", "sorna" = "sorna"))

table(scRNA_pre_filter$orig.ident)

scRNA <- scRNA_pre_filter
scRNA <- subset(scRNA, 
                subset = nFeature_RNA > 400 & nFeature_RNA < 5000 & 
                  mt_percent < 20 & 
                  HB_percent < 3 & 
                  nCount_RNA < quantile(nCount_RNA, 0.97) & 
                  nCount_RNA > 1000)

table(scRNA$orig.ident)
VlnPlot(scRNA, 
             features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"), 
             pt.size = 0.01, 
             ncol = 4, 
             group.by = "orig.ident")

theme(axis.title.x = element_text(face="bold")) + 
  scale_x_discrete(labels = c("sgrna" = "sgrna", "sorna" = "sorna"))



scRNA <- NormalizeData(scRNA) %>% # 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = T)
DimPlot(scRNA,reduction = "pca",group.by = "orig.ident")


library(harmony)

scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")

ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")

library(patchwork)

for (i in c(5,10,15, 20, 25, 30)) {
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.2)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5) + labs(title = paste0("dims: ", i)))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- FindClusters(
  object = scRNA_harmony,
  resolution = c( seq( 0.1, 1.5, 0.1) ) 
)

clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.")

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.2", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.3", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.4", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5", 
              label = T, repel = F, shuffle = T )
#
p6 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.6", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.7", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.8", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.9", 
              label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1", 
               label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.1", 
               label = T, repel = F, shuffle = T )

DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.2", 
               label = T, repel = F, shuffle = T )
DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.3", 
               label = T, repel = F, shuffle = T )
DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.4", 
               label = T, repel = F, shuffle = T )
DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.5", 
               label = T, repel = F, shuffle = T )

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.7)


scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)

umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)



library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)


markers <- c("Mki67", "Top2a","tdTomato","Sox2Gli3a","Sox2","Ascl1","Nes",# Astrocyte_lineage
             "Dcx","Dbx2","Foxp1","Plcb4","Notch1","Sox21","Ebf1","Pou3f3","Rbfox3","Rorb","Pde10a","Cux1","Tshz2","Meis2",# Astrocyte_lineage
             "Aldh1l1","Aqp4","Gfap","Col23a1","Aldoc", # Astrocyte_lineage
             "Cspg4","Pdgfra","C1ql1","Gpr17", # OPC/NG2
             "Tnr","Ust", # Oligo_P1_lineage
             "Mbp","Mobp","Ptgds","Mog","Plp1", # Oligo_P2_lineage
             "Bicc1","Dnah12","Cfap43", # Ependymal
             "Nkg7","Cd3g","Cd3e", # T_cell
             "Tyrobp","Ctss","C1qa" # Micro/Macro#F3F0F9
)

missing_genes <- markers[!markers %in% rownames(scRNA_harmony)]
print(missing_genes)
valid_markers <- markers[markers %in% rownames(scRNA_harmony)]
unique_markers <- unique(valid_markers)

dotplot <- DotPlot(
  scRNA_harmony,
  features = valid_markers,
  group.by = "seurat_clusters", 
  dot.scale = 6,         
  scale = TRUE,          
  cols = c('#008D83', '#C34A36')
) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(face = "italic")
  ) +
  labs(x = "Genes", y = "Cell Types") +
  scale_color_gradientn(
    colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'),
    values = scales::rescale(c(0, 0.5, 1))  
  )


setwd("E:\\single-cell\\data")
# 拟时序分析
library(Seurat)
#没有安装包的需要先安装这些包
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(Biobase)
library(ggsci)
library(ggpubr)
library(data.table)
library(monocle)
set.seed(12345)

load("scdata1.Rdata")


CellDimPlot3D(srt = sce, group.by = "seurat_clusters")

sce <- RunDEtest(srt = sce, group_by = "seurat_clusters", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sce, group_by = "seurat_clusters")



DEGs <- sce@tools$DEtest_celltype$Allmarkers_wilcox

DEGs <- sce@tools[["DEtest_seurat_clusters"]][["AllMarkers_wilcox"]]

sce <- NormalizeData(sce)

sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

sce <- ScaleData(sce, features = rownames(sce))

sce <- RunPCA(sce, npcs = 50)

unique(sce$celltype)
unique(scRNA_harmony$seurat_clusters)
seurat=subset(scRNA_harmony,idents = c("2","7","8","9","12","15","16","19","20"))

counts_matrix <- as(seurat[["RNA"]]$counts, "matrix")

cds <- newCellDataSet(counts_matrix,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())
library(Matrix)
library(CRsparse)

expr_matrix=seurat@assays$RNA@counts
sample_sheet<-seurat@meta.data
gene_annotation=data.frame(gene_short_name=rownames(seurat))
rownames(gene_annotation)=rownames(seurat)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))


diff_celltype <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~seurat_clusters",cores=15)
head(diff_celltype)

diff_celltype<- diff_celltype[order(diff_celltype$qval),]
ordering_genes <- row.names(diff_celltype[1:2000,])

cds <- setOrderingFilter(cds,ordering_genes = ordering_genes) 
plot_ordering_genes(cds)

cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)

p2=plot_cell_trajectory(cds, color_by = "Pseudotime")

plot_cell_trajectory(cds, color_by = "State")+
  facet_wrap(~State, nrow = 2)#
a1 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")

p3=plot_cell_trajectory(cds, color_by = "seurat_clusters")+
  facet_wrap(~seurat_clusters, nrow = 5) +
  theme(text = element_text(size = 25))  # 设置字体大小为 18

colour =c("0"="#f8766d","1"="#e68613","2"="#cd9600","3"="#aba300","4"="#7cae00","5"="#0bb70d","6"="#00be67",
  "7"="#00c19a","8"="#00bfc4","9"="#00b8e7","11"="#8494ff","10"="#00a9ff","12"="#c77cff","13"="#ed68ed","14"="#ff61cc","15"="#a58aff","16"="#d078ff","19"="#ff63b9","20"="#ff6b96"
  )
a1 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")+
  scale_color_manual(values = colour)
a1
a2 <- plot_cell_trajectory(cds, color_by = "State") 
a2
a3 <- plot_cell_trajectory(cds, color_by = "Pseudotime") 
a3
a4 <- plot_cell_trajectory(cds, color_by = "orig.ident") 
a4
#分面显示：（不同转移部位）
a5 <-plot_cell_trajectory(cds, color_by = "seurat_clusters") + facet_wrap("~seurat_clusters", nrow = 5) +
  scale_color_manual(values = colour)
a5
(a1 + a2) / (a3)
(a1 + a2) / (a3 + a4)
a6=(a1 + a2 +a3 +a4) / (a5)
a6
