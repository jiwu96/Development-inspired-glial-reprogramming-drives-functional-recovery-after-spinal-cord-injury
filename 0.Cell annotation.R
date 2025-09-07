####Data transformation and dimensionality reduction clustering############
setwd("D:/date")

library(dplyr)
library(hdf5r)
library(Seurat)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(ggplot2)

sgrna_data <- LoadH5Seurat("sgrna_filtered.h5seurat")
sorna_data <- LoadH5Seurat("sorna_filtered.h5seurat")

# Check for the presence of tdT or similar genes.
grep("tdT|tdTomato", rownames(sgrna_data), value = TRUE)

# Screen out the tdTomato-positive cells
tdT_positive_sgrna <- subset(sgrna_data, subset = tdTomato > 0)
tdT_positive_sorna <- subset(sorna_data, subset = tdTomato > 0)

sgrna <- tdT_positive_sgrna
sorna <- tdT_positive_sorna
rm(tdT_positive_sgrna,tdT_positive_sorna,sgrna_data,sorna_data)

head(sgrna@meta.data)
head(sorna@meta.data)
sgrna[["nFeature_RNA"]] <- Matrix::colSums(sgrna@assays$RNA@counts > 0)
sgrna[["nCount_RNA"]] <- Matrix::colSums(sgrna@assays$RNA@counts)
sorna[["nFeature_RNA"]] <- Matrix::colSums(sorna@assays$RNA@counts > 0)
sorna[["nCount_RNA"]] <- Matrix::colSums(sorna@assays$RNA@counts)

scRNA_pre_filter <- merge(x = sgrna, y = sorna, add.cell.ids = c("sgrna", "sorna"))

scRNA_pre_filter$orig.ident <- ifelse(grepl("sgrna_", colnames(scRNA_pre_filter)), "sgrna", "sorna")

scRNA_pre_filter[["nFeature_RNA"]] <- Matrix::colSums(scRNA_pre_filter@assays$RNA@counts > 0)
scRNA_pre_filter[["nCount_RNA"]] <- Matrix::colSums(scRNA_pre_filter@assays$RNA@counts)

scRNA_pre_filter[["mt_percent"]] <- PercentageFeatureSet(scRNA_pre_filter, pattern = "^mt-")
HB_genes_mouse <- c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-bh2", "Hbb-bh1")
HB_genes_mouse <- intersect(HB_genes_mouse, rownames(scRNA_pre_filter))
scRNA_pre_filter[["HB_percent"]] <- PercentageFeatureSet(scRNA_pre_filter, features = HB_genes_mouse)

p1 <- VlnPlot(scRNA_pre_filter, 
              features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"), 
              pt.size = 0.01, 
              ncol = 4, 
              group.by = "orig.ident")

p1 + theme(axis.title.x = element_text(face="bold")) + 
  scale_x_discrete(labels = c("sgrna" = "sgrna", "sorna" = "sorna"))
ggsave("vlnplot_front_qc.pdf", plot = p1, width = 15, height =7) 

table(scRNA_pre_filter$orig.ident)
scRNA <- scRNA_pre_filter
scRNA <- subset(scRNA, 
                subset = nFeature_RNA > 400 & nFeature_RNA < 5000 & 
                  mt_percent < 10 & 
                  HB_percent < 3 & 
                  nCount_RNA < quantile(nCount_RNA, 0.97) & 
                  nCount_RNA > 1000)

table(scRNA$orig.ident)
p <- VlnPlot(scRNA, 
             features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"), 
             pt.size = 0.01, 
             ncol = 4, 
             group.by = "orig.ident")

# 修改x轴标签
p + theme(axis.title.x = element_text(face="bold")) + 
  scale_x_discrete(labels = c("sgrna" = "sgrna", "sorna" = "sorna"))
ggsave("vlnplot_after_qc.pdf", plot = p, width = 15, height =7) 

rm(scRNA_pre_filter,sorna,sgrna,p,p1)

####Data normalization, screening of highly variable genes and PCA dimensionality reduction####

scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
a=DimPlot(scRNA,reduction = "pca",group.by = "orig.ident")
a

top15 <- head(VariableFeatures(scRNA), 15) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )

####Cell cycle score####
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))

s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
scRNA=CellCycleScoring(object = scRNA, 
                       s.features = s_genes, 
                       g2m.features = g2m_genes, 
                       set.ident = TRUE)
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
p4=VlnPlot(scRNA, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0.1)
p4
p5=scRNA@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()

ggsave(file = "s.G2M,score.pdf",plot = p4,he = 5,wi = 6 )
ggsave(file = "s.G2M,score2.pdf",plot = p5,he = 5,wi = 6 )

library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
b

pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated

ggsave(file = "pca_harmony.pdf",plot = pca_harmony_integrated,he = 5,wi = 6 )

ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")

library(patchwork)

####debugging dims####
for (i in c(5,10,15, 20, 25, 30)) {
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.2)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5) + labs(title = paste0("dims: ", i)))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_5 + plot_10 + plot_15 + plot_20 + plot_25 + plot_30 +plot_layout( ncol = 3 )
ggsave( p, filename = "as-dims.5-30.pdf", width = 15, height = 8 )
print(plot_10)
ggsave( plot_10, filename = "as-dims.10.pdf", width = 5, height = 5 )

####debugging resolution####
library(clustree)

set.seed(123)
scRNA_harmony <- scRNA_harmony_filtered
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- FindClusters(
  object = scRNA_harmony,
  resolution = c( seq( 0.1, 1.5, 0.1) ) 
)
clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.")
ggsave( filename = "clustree.pdf", width = 12, height = 9 )

p1 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", 
              label = T, repel = F, shuffle = T )
#
p2 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.2", 
              label = T, repel = F, shuffle = T )
#
p3 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.3", 
              label = T, repel = F, shuffle = T )
#
p4 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.4", 
              label = T, repel = F, shuffle = T )
#
p5 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5", 
              label = T, repel = F, shuffle = T )
#
p6 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.6", 
              label = T, repel = F, shuffle = T )
#
p7 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.7", 
              label = T, repel = F, shuffle = T )
#
p8 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.8", 
              label = T, repel = F, shuffle = T )
#
p9 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.9", 
              label = T, repel = F, shuffle = T )
#
p10 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1", 
               label = T, repel = F, shuffle = T )
p11 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.1", 
               label = T, repel = F, shuffle = T )
p12 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.2", 
               label = T, repel = F, shuffle = T )
p13 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.3", 
               label = T, repel = F, shuffle = T )
p14 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.4", 
               label = T, repel = F, shuffle = T )
p15 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1.5", 
               label = T, repel = F, shuffle = T )

p = p1 + p2 + p3 + p4 + p5 + p6 + p7+ p8+ p9+ p10 +p11 +p12 +p13 +p14 +p15 +plot_layout( ncol = 5 )
ggsave( p, filename = "as-RNA_snn_res.0.1-1.5.pdf", width = 25, height = 15 )

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.8)
table(scRNA_harmony@meta.data$seurat_clusters)

##umap/tsne
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)


umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)

umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)

umap_tsne_integrated

ggsave("as-10-dims_umap_tsne_integrated1.pdf",umap_tsne_integrated,wi=15,he=15)

save(scRNA_harmony,file = "as_dims.Rdata")

rm(list = ls())
load("as_dims.Rdata")


########scRNA_harmony 
####FindAllMarkers
View(scRNA_harmony@meta.data)
table(scRNA_harmony@meta.data$orig.ident)

# 然后运行 FindAllMarkers
markers <- FindAllMarkers(
  object = scRNA_harmony, 
  test.use = "wilcox", 
  only.pos = TRUE,
  logfc.threshold = 0.25
)

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20,"cluster_top20.csv",row.names = T)

markers <- c(
  'Clec4c', 'Il3ra', 'Nrp1', 'Cd1c', 'Cst3', 'Fcer1a','Ptprc', 'Cd14', 
  'Mki67', 'Stmn1', 'Top2a', 'S100a8', 'S100a9', 'Vwf', 'Cldn5', 'Cdh5', 
  'Col1a1', 'Dcn', 'C1qa', 'Tubb3', 'Fcgr1', 'Csf1r', 'Nos2', 'Il12b', 
  'Cxcl9', 'Cxcl10', 'Cd80', 'Cd86', 'Cd206', 'Mrc1', 'Il10', 'Arg1', 
  'Cd209a', 'Flt3', 'Xcr1', 'Itgae', 'Itgam', 'Ccr9', 'Zbtb46', 'Nudt17', 
  'Itgax', 'Cd3e', 'Lat', 'Bcl11b', 'Gata3', 'Il1rl1', 'Hs3st1', 'Rnf128', 
  'Il2ra', 'Acta2', 'Adgre1', 'Ly6c2', 'Ccr2', 'Chil3', 'Ace', 'Eno3', 
  'H2-Ab1', 'Sparc', 'Aldh1l1', 'Gfap', 'Slc1a3', 'Cspg4', 'Pdgfra', 
  'Olig1', 'Olig2', 'Mog', 'Mbp', 'Ascl1', 'Sox2', 'Sox10', 'Ngn2', 
  'Dcx', 'Rbfox3', 'Map2', 'Syp', 'Zfhx3', 'Adamts2', 'Pitx2', 'Gpc3', 
  'Sema5b', 'Ccbe1', 'Tfap2b', 'Gad2'
)
DotPlot(scRNA_harmony, features = markers)
DotPlot(scRNA_harmony, features = markers) + coord_flip() +  
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +  
  theme(axis.text.x = element_text(angle = 90)) +  
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

ggsave(filename="GS-mark1.pdf", width = 350, height = 300, units = "mm")


#####16cluster_oligodendroglia
celltype=data.frame(ClusterID=0:13,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='Astro-lineage cells'
celltype[celltype$ClusterID %in% c(1),2]='OPCING2 lineage'
celltype[celltype$ClusterID %in% c(2),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(3),2]='Astro-lineage cells'
celltype[celltype$ClusterID %in% c(4),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(5),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(6),2]='Astro-lineage cells'
celltype[celltype$ClusterID %in% c(7),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(8),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(9),2]='Micro/Macro'
celltype[celltype$ClusterID %in% c(10),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(11),2]='Oligodendrocytes'
celltype[celltype$ClusterID %in% c(12),2]='Ependymal cell'
celltype[celltype$ClusterID %in% c(13),2]='T_cell'

celltype
table(celltype$celltype)
sce.in=scRNA_harmony
sce.in@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.in@meta.data$celltype)


table(sce.in@meta.data$celltype,sce.in@meta.data$seurat_clusters)
sce=sce.in
p.dim.cell=DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T,pt.size = 1) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="DimPlot_tsne_celltype.pdf",width=9, height=7)
p.dim.cell=DimPlot(sce, reduction = "umap", group.by = "celltype",label = F,pt.size = 1,repel = T) 
p.dim.cell
library(ggplot2)
ggsave(plot=p.dim.cell,filename="DimPlot_umap_celltype.pdf",width=9, height=7)


VlnPlot(sce, features = "Dcx",group.by = "celltype",pt.size = 1)
VlnPlot(sce, features = "SOX9",group.by = "celltype",pt.size = 1)
DefaultAssay(sce) <- "RNA"
#
save(sce,file='scdata.Rdata')
