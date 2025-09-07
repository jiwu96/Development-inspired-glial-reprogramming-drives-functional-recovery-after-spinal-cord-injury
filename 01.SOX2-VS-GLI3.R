####SOX2-VS-GLI3 fig1

setwd("D:/date")
load("scdata.Rdata")

scRNA_harmony=sce

scRNA_harmony@meta.data$seurat_clusters <- as.character(scRNA_harmony@meta.data$seurat_clusters)

scRNA_harmony@meta.data$seurat_clusters <- scRNA_harmony@meta.data$seurat_clusters

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("6","11","4","10","15","9","17","0","13")
] <- "Neurogeneration pool"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("3","12")
] <- "OPC/NG2_starter"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("7","8")
] <- "LVinfectedOligo_P"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("16","5","2")
] <- "Oligo_P"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("19","1","21")
] <- "Oligodendrocyte"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("14")
] <- "Micro/Macro"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("18")
] <- "Ependymal_cell"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("20")
] <- "Tcell"

scRNA_harmony@meta.data$seurat_clusters[
  scRNA_harmony@meta.data$seurat_clusters %in% c("22")
] <- "NA"

table(scRNA_harmony$seurat_clusters)

head(scRNA_harmony@meta.data$seurat_clusters)

library(Seurat)
library(SeuratData)
library(clusterProfiler)
library(org.Mm.eg.db)

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(scRNAtoolVis)
library(devtools)
library(reshape2)
library(ggplot2)
library(tidydr)
library(RColorBrewer)
scRNA=scRNA_harmony
scRNA=scRNA_harmony
DimPlot=DimPlot(scRNA, reduction = "umap", group.by = "orig.ident",label = F,pt.size = 1,repel = T)
DimPlot
DimPlot=DimPlot(scRNA, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T)
DimPlot

ggsave(filename = "cluster-DimPlot1.pdf", plot = DimPlot, width = 350, height = 250, units = "mm")

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

dotplot <- DotPlot(scRNA, features = markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))
dotplot


####按照自己想的细胞顺序进行排序
# Ensure the correct order of clusters
scRNA_harmony$seurat_clusters <- factor(scRNA_harmony$seurat_clusters,
                                        levels = c("Neurogeneration pool",
                                                   "OPC/NG2_starter",
                                                   "LVinfectedOligo_P",
                                                   "Oligo_P",
                                                   "Oligodendrocyte",
                                                   "Ependymal_cell",
                                                   "Tcell",
                                                   "Micro/Macro",
                                                   "NA"))

# Generate the DotPlot with the specified order of clusters
dotplot <- DotPlot(scRNA_harmony, features = markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

# Plot the result
dotplot
ggsave(filename = "25-01-06-9cluster-marker.pdf", plot = dotplot, width = 115, height = 250, units = "mm")



####delete NA####
selected_clusters <- c(
  "Neurogeneration pool",
  "OPC/NG2_starter",
  "LVinfectedOligo_P",
  "Oligo_P",
  "Oligodendrocyte",
  "Ependymal_cell",
  "Micro/Macro",
  "Tcell"
)

scRNA1 <- subset(scRNA, seurat_clusters %in% selected_clusters)

#####Multiple sets of volcano maps########
scRNA=scRNA1
unique(Idents(scRNA))
Idents(scRNA) <- "new_cluster"
Idents(scRNA) <- "seurat_clusters"
table(scRNA$new_cluster)
table(scRNA$orig.ident)

# 加载包
library(Seurat)
library(tidyverse)
library(ggrepel)
library(scRNAtoolVis)

# 1. 提取三群
selected_groups <- subset(scRNA, idents = c("Neurogeneration pool",
                                            "LVinfectedOligo_P"))

# 确保 orig.ident 是分组变量new_cluster
Idents(selected_groups) <- selected_groups$new_cluster

# 2. 对每个 new_cluster 分别做 sgrna vs sorna 差异分析
deg_list <- list()

for (cl in c("Neurogeneration pool", "OPC/NG2_starter", "LVinfectedOligo_P")) {
  tmp <- subset(selected_groups, idents = c("sgrna", "sorna"))
  tmp <- subset(tmp, subset = new_cluster == cl)
  Idents(tmp) <- tmp$orig.ident
  
  # 差异分析
  deg_tmp <- FindMarkers(tmp,
                         ident.1 = "sgrna",
                         ident.2 = "sorna",
                         min.pct = 0.25,
                         logfc.threshold = 0.25)
  
  deg_tmp$cluster <- cl
  deg_tmp$gene <- rownames(deg_tmp)
  deg_list[[cl]] <- deg_tmp
}

# 合并差异结果
deg_all <- bind_rows(deg_list)

# 3. 保存差异分析结果
write.table(deg_all, file = "多组火山图.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)
colnames(deg_all)
head(deg_all)
table(deg_all$cluster)

dd <- deg_all %>%
  transmute(
    gene,
    cluster = as.character(cluster),
    logFC = avg_log2FC,
    pvalue = p_val
  ) %>%
  filter(is.finite(logFC), is.finite(pvalue)) %>%
  mutate(cluster = factor(cluster, levels = c("Neurogeneration pool",
                                              "OPC/NG2_starter",
                                              "LVinfectedOligo_P")))
dd <- dd %>%
  mutate(regulate = case_when(
    pvalue < 0.05 & logFC >= 0.585 ~ "Up",
    pvalue < 0.05 & logFC <= -0.585 ~ "Down",
    TRUE ~ "Not"
  ))

tile.cols <- c("#8DD3C7", "#FFFFB3", "#80B1D3") # 与因子顺序一一对应
stopifnot(nlevels(dd$cluster) == length(tile.cols))



# 改成 jjVolcano 能识别的列名
deg_all <- deg_all %>%
  rename(avg_log2FC = avg_log2FC,  # 如果已经是这个名字就不用改
         p_val = p_val) %>%
  mutate(cluster = as.character(cluster))

# 过滤检查
sum(deg_all$p_val < 0.05 & abs(deg_all$avg_log2FC) >= 0.585)


# 4. 生成多组火山图
colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

pdf("多组火山图3组.pdf", width = 12, height = 6)
jjVolcano(
  diffData = deg_all,
  tile.col = colors[1:length(unique(deg_all$cluster))], # 每组底色
  aesCol = c('purple','orange'), # 下调/上调颜色
  pSize = 0.8,                    # 点大小
  pvalue.cutoff = 0.05,           # P值阈值
  log2FC.cutoff = 1,              # log2FC阈值
  topGeneN = 3,                   # 每组标注的top基因数
  col.type = "upDown"             # 上下调配色
)
dev.off()


####Cell ratio chart#####

colour=c("#53b400","#f98b84","#a58aff","#00c094",
         "#fb61d7", "#c49a00","#00b6eb")



colour <- c("0" = "#faaca7",
            "10" = "#e1c066",
            "11" = "#ffa0e0",
            "13" = "#ddb0ff",
            "15" = "#66cbff",
            "17" = "#66d8db",
            "4" = "#66d8a0",
            "9" = "#b0ce66" )

#p1 <- plot_cell_trajectory(cds, color_by = "seurat_clusters") +

scRNA=scRNA_harmony_filtered

Idents(scRNA) <- scRNA$new_cluster

prop.table(table(Idents(scRNA)))

table(Idents(scRNA), scRNA$orig.ident)


Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("seurat_clusters", "sample", "ratio")


Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("seurat_clusters", "sample", "ratio")


ggplot(Cellratio, aes(x = sample, y = ratio, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "fill", colour = NA) +  # 去除线条
  scale_fill_manual(values = colour) +  # 使用自定义颜色
  labs(title = "Cell Type Proportions by Sample",
       x = "Sample",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


table(scRNA$orig.ident)
prop.table(table(Idents(scRNA)))  
table(Idents(scRNA), scRNA$orig.ident)  
Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)  
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("seurat_clusters","sample","ratio")



colourCount = length(unique(Cellratio$seurat_clusters)) 
ggplot(Cellratio) + 
  geom_bar(aes(x = sample, y = ratio, fill = seurat_clusters), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +  # 添加条形图层，根据细胞类型填充颜色
  theme_classic() + 
  labs(x = 'Sample', y = 'Ratio') +  
  coord_flip() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))+   # 设置面板边框属性
  scale_fill_manual(values = colour)  

ggplot(Cellratio) +
  geom_bar(aes(x =sample, y= ratio, fill = seurat_clusters),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = colour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


#####单个gene在umap图上的表达#######

scRNA1 <- subset(scRNA, seurat_clusters %in% selected_clusters)

scRNA=scRNA1
scRNAsub=scRNA
unique(Idents(scRNA))
Idents(scRNA) <- "seurat_clusters"


plot=FeaturePlot(scRNAsub, features = c('tdTomato',"Sox2Gli3a","Ascl1","Top2a","Gfap","Dcx","Rbfox3","Dbx2","C1qa","Ust"), pt.size = 1, order = TRUE)
plot


plot=FeaturePlot(scRNAsub, features = c('Sox2','tdTomato'), pt.size = 1, order = TRUE)
plot

plot=FeaturePlot(scRNAsub, features = c('Cdkn1a','Cdkn1b','Cdk4','Cdk6','Pcna','Mcm2','Mcm5','Cdk1','Ccnb1','Ccna2','Tubb','Ccnb1','Plk1'), pt.size = 1, order = TRUE)
plot
ggsave(filename = "Sox2-tdt-Sox2Gli3a-marker1.pdf", plot = plot, width = 455, height = 400, units = "mm")


library(Nebulosa)
library(ggrastr)

density_plot <- ggrastr::rasterize(
  Nebulosa::plot_density(scRNAsub, c('tdTomato',"Sox2Gli3a","Ascl1","Top2a","Gfap","Dcx","Rbfox3","Dbx2","C1qa","Ust"), size = 0.2),
  dpi = 300
)

# 显示绘图
print(density_plot)

# 保存图片为高分辨率的 PDF 文件
ggsave("Sox2Gli3a.pdf", plot = density_plot, device = "pdf", width = 10, height = 8)
FeaturePlot

# 如果你有多个降维结果，请确保使用的是正确的降维类型
DimPlot(scRNAsub, reduction = "umap")

DimPlot(scRNAsub, label=TRUE)
FeaturePlot(scRNAsub,features = 'Sox2Gli3a',pt.size = 1,order = T)
p2=FeaturePlot(scRNAsub, features = 'Ust', pt.size = 1, order = TRUE) +
  scale_color_gradientn(colors = c("blue", "yellow", "red")) +  # Custom color gradient for expression levels
  theme_minimal()  # Use a clean theme for better visualization
p2
ggsave(filename = "Ust.pdf", plot = p2, width = 300, height = 250, units = "mm")

dev.off()

######fig2-J图########
load("Neurogeneration-NG2-LV-infected_Oligo_P_cluster.Rdata")
markers <- c("Mafa","Sp8","Pou6f2","Foxp2")#,"","",
#####12-18 神经元marker###### "Cdk4","Pcna","Eomes"#,"",
markers <- c("Gfap","Aldoc","Ebf1","Pou3f3",
             "Meis2","Rtn1","Map2","Erbb4",
             "Pde10a","Tshz2","Gria2","Rorb","Cux1","Thsd7a","Prox1","Sox1",
             "Notch1","Zfpm2","Cmip","Sox21","Foxp1","Foxp2","Plcb4","Rbfox3","Sox5",
             "Hopx","Egfr","Olig1","Olig2","Nkx2-2","Nkx6-2","Dbx2","Ascl1",
             "Sox2","Cdk4","Pcna","Mki67","Top2a","Sox2Gli3a","tdTomato")
#"Hoxa1","Hoxa2","Hoxa3","Hoxa4","Hoxa5","Hoxa6","Hoxa7","Hoxa9","Hoxa10","Hoxa11","Hoxa13","Hoxb1","Hoxb2","Hoxb3","Hoxb4","Hoxb5","Hoxb6","Hoxb7","Hoxb8","Hoxb9","Hoxb13","Hoxc4","Hoxc5","Hoxc6","Hoxc8","Hoxc9","Hoxc10","Hoxc11","Hoxc12","Hoxc13","Hoxd1","Hoxd3","Hoxd4","Hoxd8","Hoxd9","Hoxd10","Hoxd11","Hoxd12","Hoxd13")

# Check if scRNA is a Seurat object
scRNA_harmony=scRNA
# 检查当前的 cluster 标签
unique(scRNA$seurat_clusters)

# 如果有 NA 值，先移除它们
#scRNA <- subset(scRNA, subset = !is.na(seurat_clusters))

# 再次检查
unique(scRNA$seurat_clusters)

# 定义你想要的顺序
desired_order <- c("6", "9",  "15", "4", "10", "17", "0", "13", "11","12","3","8","7")

# 检查所有指定的 levels 是否存在于数据中
missing_levels <- setdiff(desired_order, unique(scRNA$seurat_clusters))
if (length(missing_levels) > 0) {
  warning("以下 cluster 不存在于数据中: ", paste(missing_levels, collapse = ", "))
}

# 重新定义因子顺序
scRNA$seurat_clusters <- factor(scRNA$seurat_clusters, levels = desired_order)

# 验证顺序是否正确
levels(scRNA$seurat_clusters)

# 检查 markers 是否存在于数据中
valid_markers <- markers[markers %in% rownames(scRNA)]

# 生成 DotPlot
dotplot <- DotPlot(scRNA, features = valid_markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36')
  )
dotplot
# 显示图形
print(dotplot)

head(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$orig.ident)
# 合并 cluster 和 orig.ident 信息
scRNA$cluster_group <- paste0(scRNA$seurat_clusters, "_", scRNA$orig.ident)

# 重新设置顺序（保证 cluster 顺序不乱）
scRNA$cluster_group <- factor(
  scRNA$cluster_group,
  levels = unlist(lapply(desired_order, function(cl) paste0(cl, "_", c("sgrna", "sorna"))))
)

# 生成 DotPlot
dotplot <- DotPlot(scRNA, features = valid_markers, group.by = "cluster_group") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36')
  )
dotplot


ggsave(filename = "8-28GLI3-all-cluster-marker.pdf", plot = dotplot, width = 220, height = 230, units = "mm")

dotplot <- DotPlot(scRNA, features = valid_markers,
                   group.by = "seurat_clusters", split.by = "orig.ident") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3))
dotplot






scRNA$seurat_clusters <- factor(scRNA$seurat_clusters, levels = c("6", "9", "12", "15", "4", "10", "17", "0", "13", "11"))

# Generate DotPlot with valid markers
dotplot <- DotPlot(scRNA, features = markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

print(dotplot)


ggsave(filename = "8-15 GLI3-all-cluster-marker_updated2.pdf", plot = dotplot, width = 125, height = 240, units = "mm")



####KEGG-go######

rm(list=ls())
library(Seurat)
library(SeuratData)
library(clusterProfiler)
library(org.Mm.eg.db)
sce=selected_groups
sce=scRNA
head(sce@meta.data)

Idents(sce) <- "new_cluster"

scRNA1 <- subset(sce, idents = c("Neurogeneration pool", "LVinfectedOligo_P"))

# 差异分析
library(Seurat)
diff_markers <- FindMarkers(scRNA1, ident.1 = "Neurogeneration pool", ident.2 = "LVinfectedOligo_P")

head(diff_markers)

library(clusterProfiler)
library(org.Mm.eg.db)

significant_genes <- rownames(diff_markers[diff_markers$p_val_adj < 0.05, ])

gene_entrez <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

options(timeout = 300)  # 将超时时间设置为 300 秒
kegg_results <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "mmu")
kegg_results <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "mmu")

head(kegg_results)


#####提取 sgrna 和 sorna#####
sce_group <- subset(scRNA1, orig.ident %in% c("sgrna", "sorna"))


group_markers <- FindMarkers(sce_group, ident.1 = "sgrna", ident.2 = "sorna", group.by = "orig.ident")

head(group_markers)

group_significant_genes <- rownames(group_markers[group_markers$p_val_adj < 0.05, ])

group_gene_entrez <- bitr(group_significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

group_kegg_results <- enrichKEGG(gene = group_gene_entrez$ENTREZID, organism = "mmu")

head(group_kegg_results)

library(ggplot2)

# Neurogeneration pool 和 OPC/NG2_starter 的 KEGG 结果
barplot(kegg_results, showCategory = 10, title = "KEGG Pathways (Neurogeneration vs OPC/NG2)")

# sgrna 和 sorna 的 KEGG 结果
barplot(group_kegg_results, showCategory = 10, title = "KEGG Pathways (sgrna vs sorna)")




# 得到每个细胞亚群的高表达基因
obj.markers=significant_genes
obj.markers <- FindAllMarkers(sce, only.pos = TRUE)
head(obj.markers)

Symbol <- mapIds(get("org.Mm.eg.db"), keys = obj.markers$gene, keytype = "SYMBOL", column="ENTREZID")
ids <- bitr(obj.markers$gene,"SYMBOL","ENTREZID", "org.Mm.eg.db")

# 合并ENTREZID到obj.markers中
data <- merge(obj.markers, ids, by.x="gene", by.y="SYMBOL")

head(data)


gcSample <- split(data$ENTREZID, data$cluster)
gcSample
# Perform GO enrichment analysis
xx_go <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Mm.eg.db", ont="BP", pvalueCutoff=1, qvalueCutoff=1)

# Perform KEGG enrichment analysis
xx_kegg <- compareCluster(gcSample, fun="enrichKEGG", organism="mmu", pvalueCutoff=1, qvalueCutoff=1)



res <- xx_go@compareClusterResult

## 将富集结果中的 ENTREZID 重新转为 SYMBOL
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}

head(res)

library(dplyr)
## 通路筛选
enrich <- res %>%
  group_by(Cluster) %>%
  top_n(n = 5, wt = -pvalue) %>%
  filter(Cluster %in% c("Astrocyte-lineage", "OPC/NG2-lineage"))

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
colnames(dt)

library(ggplot2)
library(ggridges)
# 先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)

p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "OPC/NG2-lineage         Astrocyte-lineage", title = "GO BP enrichment") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, color=rep(c('#6bb9d2', '#d55640'),each=5)) + # hjust = 0,左对齐
  theme_classic() +
  mytheme +
  NoLegend()

p



res <- xx_kegg@compareClusterResult

## 将富集结果中的 ENTREZID 重新转为 SYMBOL
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}

head(res)

library(dplyr)
## 通路筛选
enrich <- res %>%
  group_by(Cluster) %>%
  top_n(n = 5, wt = -pvalue) %>%
  filter(Cluster %in% c("Astrocyte-lineage", "OPC/NG2-lineage"))

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
colnames(dt)

library(ggplot2)
library(ggridges)
# 先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)

p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "OPC/NG2-lineage         Astrocyte-lineage", title = "KEGG Pathway enrichment") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, color=rep(c('#6bb9d2', '#d55640'),each=5)) + # hjust = 0,左对齐
  theme_classic() +
  mytheme +
  NoLegend()

p

# 保存，这里的保存宽和高进行了调整，可以使得结果比较美观
ggsave(filename = "kegg-test.png", width = 5.2, height = 5.1, plot = p)

####KEGG-go美化######

library(dplyr)
markers <- sce.markers |> group_by(cluster) |>
  filter(p_val_adj < 0.001) |>
  ungroup()

library(clusterProfiler)
gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Mm.eg.db')
markers <- full_join(markers, gid, by=c('gene' = 'SYMBOL'))
xx_kegg = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG')

dotplot(xx_kegg, label_format=40) + theme(axis.text.x = element_text(angle=45, hjust=1))



library(COSG)
marker_cosg <- cosg(
  sce,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)

head(marker_cosg[[1]])

y = compareCluster(marker_cosg[[1]], fun='enrichGO',
                   OrgDb = 'org.Mm.eg.db', keyType = 'SYMBOL', ont="MF")

dotplot(xx_go, label_format=60) + theme(axis.text.x = element_text(angle=45, hjust=1))



######Neurogeneration VS OPC#########
# 确定分组
Idents(sce) <- "seurat_clusters"

# 提取差异基因
deg_seurat_clusters <- FindMarkers(
  object = sce,
  ident.1 = "Neurogeneration pool",
  ident.2 = "OPC/NG2_starter",
  logfc.threshold = 0.25, # 可调整阈值
  min.pct = 0.1          # 可调整阈值
)
head(deg_seurat_clusters)

# 筛选 adj.p.val < 0.05 且 logFC > 0 的基因
sig_genes <- rownames(deg_seurat_clusters[deg_seurat_clusters$p_val_adj < 0.05 & abs(deg_seurat_clusters$avg_log2FC) > 0.25, ])

# 基因名转化为 ENTREZ ID
gene_entrez <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

###GO富集分析
go_results <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",        # 可选 BP, CC, MF
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)
head(go_results)

##kegg富集分析
options(clusterProfiler.download.method = "wininet")
kegg_results <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
head(kegg_results)

barplot(go_results, showCategory = 10, title = "GO Enrichment Analysis")

dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")
#####################


#####1.kegg-go富集分析#####
library(ggplot2)
de_markers=diff_markers
de_markers=deg_seurat_clusters
# Add a column to identify significant genes based on thresholds (e.g., p-value < 0.05, log2FC > 0.5)
de_markers$gene <- rownames(de_markers)
de_markers$significance <- ifelse(de_markers$p_val_adj < 0.05 & abs(de_markers$avg_log2FC) > 1,
                                  ifelse(de_markers$avg_log2FC > 1, "Upregulated", "Downregulated"), "Not Significant")

# Volcano plot
# 绘制火山图
p <- ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_minimal() +
  # 设置自定义颜色
  scale_color_manual(values = c("Upregulated" = "#C34A36", "Downregulated" = "#6bb9d2", "Not Significant" = "gray")) +
  # 设置标题
  ggtitle("Differential Gene Expression - Neurogeneration pool VS Neural precursor cells") +
  # 显示 x 和 y 轴
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  # 去掉中间网格，仅保留边框
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )
p
# 保存为 PDF
ggsave("Neurogeneration_pool_VS_Neural precursor cells-vol.pdf", plot = p, width = 8, height = 6)


# Extract upregulated and downregulated genes
upregulated_genes <- rownames(subset(de_markers, avg_log2FC > 1 & p_val_adj < 0.05))
downregulated_genes <- rownames(subset(de_markers, avg_log2FC < 1 & p_val_adj < 0.05))

# All differentially expressed genes
all_genes <- rownames(subset(de_markers, p_val_adj < 0.05))
####go-kegg富集######

library(clusterProfiler)
library(org.Mm.eg.db)  # For human gene annotations, or org.Mm.eg.db for mouse
library(enrichplot)    # For enrichment visualization

# GO enrichment for upregulated genes
go_up <- enrichGO(gene = upregulated_genes,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

# GO enrichment for downregulated genes
go_down <- enrichGO(gene = downregulated_genes,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)

# GO enrichment for all DE genes
go_all <- enrichGO(gene = all_genes,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

##图片合并并美化

library(dplyr)
library(clusterProfiler)
library(ggplot2)

# 提取 GO 结果
go_up_df <- as.data.frame(go_up) %>%
  mutate(Group = "Upregulated")
go_down_df <- as.data.frame(go_down) %>%
  mutate(Group = "Downregulated")

###go_up可视化####
# 筛选每组前 top N（例如前 10 个 term）
top_n <- 10
go_up <- go_up %>%
  group_by(Group) %>%
  slice_min(order_by = p.adjust, n = top_n) %>%
  ungroup()

# 添加颜色信息
go_up <- go_up %>%
  mutate(Color = ifelse(Group == "Upregulated", "#C34A36"))

# 绘图
p <- ggplot(go_up, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = Color)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  facet_wrap(~Group, scales = "free_y") +  # 按组分面
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 移除 x 轴网格
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"), # 仅保留 y 轴网格
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_identity() +  # 使用自定义颜色
  labs(x = "Gene Ratio", y = NULL, title = "GO Enrichment for Up- and Downregulated Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

# 显示图表
print(p)



###go_down可视化####
# 筛选每组前 top N（例如前 10 个 term）
top_n <- 10
go_down <- go_down %>%
  group_by(Group) %>%
  slice_min(order_by = p.adjust, n = top_n) %>%
  ungroup()

# 添加颜色信息
go_down <- go_down %>%
  mutate(Color = ifelse(Group == "Upregulated", "#C34A36", "#6bb9d2"))

# 绘图
p <- ggplot(go_down, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = Color)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  facet_wrap(~Group, scales = "free_y") +  # 按组分面
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 移除 x 轴网格
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"), # 仅保留 y 轴网格
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_identity() +  # 使用自定义颜色
  labs(x = "Gene Ratio", y = NULL, title = "GO Enrichment for Up- and Downregulated Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

# 显示图表
print(p)


# 合并数据
go_combined <- bind_rows(go_up_df, go_down_df)

# 筛选每组前 top N（例如前 10 个 term）
top_n <- 10
go_combined <- go_combined %>%
  group_by(Group) %>%
  slice_min(order_by = p.adjust, n = top_n) %>%
  ungroup()

# 添加颜色信息
go_combined <- go_combined %>%
  mutate(Color = ifelse(Group == "Upregulated", "#C34A36", "#6bb9d2"))

# 绘图
p <- ggplot(go_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = Color)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  facet_wrap(~Group, scales = "free_y") +  # 按组分面
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 移除 x 轴网格
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"), # 仅保留 y 轴网格
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_identity() +  # 使用自定义颜色
  labs(x = "Gene Ratio", y = NULL, title = "GO Enrichment for Up- and Downregulated Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

# 显示图表
print(p)




####Suggested chronological analysis#####
library(ClusterGVis)
library(Seurat)
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
setwd("D:\\date")
load("cluster细胞注释.Rdata")
sce=scRNA_harmony
Idents(sce)=sce$ seurat_clusters
DimPlot(sce)

####创建monocle对象####
unique(sce$ seurat_clusters)
seurat=subset(sce,idents = c("19","1","21","16","5","2","12","3","7","8","6","4","10","15","9","17","0","13","11"))
unique(seurat$ seurat_clusters)
table(seurat$ seurat_clusters)

expr_matrix=seurat@assays$RNA@counts
sample_sheet<-seurat@meta.data
gene_annotation=data.frame(gene_short_name=rownames(seurat))
rownames(gene_annotation)=rownames(seurat)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())
cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_celltype <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~seurat_clusters",cores=12)#差异分析
head(diff_celltype)#查看前几行数据
write.csv(diff_celltype,'all-Neuro_degForCellOrdering.csv')


diff_celltype<- diff_celltype[order(diff_celltype$qval),]
ordering_genes <- row.names(diff_celltype[1:1000,])

cds <- setOrderingFilter(cds,ordering_genes = ordering_genes) 
plot_ordering_genes(cds)
ggsave(filename = 'all-Neuromonocle2_ordering_gene.pdf')

#
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)

p1=plot_cell_trajectory(cds, color_by = "State") +
  theme(text = element_text(size = 18))  
p1
ggsave(p1,filename = 'all-Neuromonocle2_state_trajectory.pdf', width = 12, height = 9)
p1=plot_cell_trajectory(cds, color_by = "seurat_clusters") +
  theme(text = element_text(size = 18)) 

p1
ggsave(p1,filename = 'sgrna-Neuromonocle2_seurat_clusters_trajectory.pdf', width = 12, height = 9)


#4.按发育时间分类画轨迹图，可以与state的图对应上
p2=plot_cell_trajectory(cds, color_by = "Pseudotime")+
  theme(text = element_text(size = 18))  
p2
ggsave(p2,filename = 'allNeuro-monocle2_Pseudotime.pdf', width = 12, height = 9)





library(ClusterGVis)

clusterAndVisualize <- function(cds, num_clusters = 12, gene_num = 2000, enrich_db = "org.Mm.eg.db",
                                organism = "mmu", pvalue_cutoff = 0.05, topn = 5, seed = 123456,
                                output_pdf = "拟时序生物学变化1.pdf",
                                pdf_height = 20, pdf_width = 12, output_table = "拟时序clustered_genes.csv") {
  # 差异基因测试
  diff_test <- differentialGeneTest(cds, cores = 16, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_genes <- diff_test %>% arrange(qval) %>% head(gene_num) %>% dplyr::select(gene_short_name)
  
  df <- plot_pseudotime_heatmap2(cds[diff_genes[,1], ],
                                 num_clusters = num_clusters,
                                 cores = 1,
                                 use_gene_short_name = TRUE,
                                 show_rownames = TRUE)
  
  # 富集分析
  enrich <- enrichCluster(object = df,
                          OrgDb = enrich_db,
                          type = "BP",
                          organism = organism,
                          pvalueCutoff = pvalue_cutoff,
                          topn = topn,
                          seed = seed)
  
  # 创建PDF输出
  pdf(file = output_pdf, height = pdf_height, width = pdf_width, onefile = FALSE)
  
  # 可视化聚类结果和富集分析结果
  visCluster(object = df,
             plot.type = "both",
             column_names_rot = 20,
             show_row_dend = FALSE,
             markGenes = sample(df$wide.res$gene, gene_num, replace = FALSE),
             markGenes.side = "right",
             annoTerm.data = enrich,
             go.col = rep(jjAnno::useMyCol("calm", n = num_clusters), each = 5),
             line.side = "left")
  
  # 关闭PDF设备
  dev.off()
}

# 使用修改的人类数据进行聚类和可视化
clusterAndVisualize(cds, num_clusters = 12, gene_num = 2000, enrich_db = "org.Mm.eg.db", organism = "mmu")


#####绘制 3D 图######
library(Seurat)
library(dplyr)
library(plotly)

af=sce
DimPlot(af, pt.size = 0.2, group.by = "seurat_clusters", reduction = "umap")
# 2. 3D scatter plot

af <- RunUMAP(af, dims = 1:10, reduction = "pca", n.components = 3)
af <- RunTSNE(af, dims = 1:10, reduction = "pca", n.components = 3)
umap <- af@reductions$umap@cell.embeddings %>% as.data.frame()
umap$seurat_clusters <- af$seurat_clusters %>% as.character()

# 提取 UMAP 嵌入
umap <- af@reductions$umap@cell.embeddings %>% as.data.frame()
umap$seurat_clusters <- af$seurat_clusters %>% as.character()


cluster_colors <- c("9"="#fe6a9a","10"= "#d59100","17"="#00c08e","0"= "#f8766c","3"="#00a7fe",
                    "4"="#7e95fe","8"="#fe61be","1"= "#e9842c","7"="#f762df","2"="#00bcd3",
                    "16"="#00bc60", "11"="#bb9d00","15"="#00b813","13"="#6eaf00","6"="#e16df7",
                    "5"="#bb81fe",  "12"="#9ca700", "21"="#00b5ed","19" = "#00c0b4"
)

# 为 UMAP 数据添加颜色
umap <- umap %>%
  mutate(
    Color = cluster_colors[as.character(seurat_clusters)]
  )

# 绘制 3D 图
plot_3d <- plot_ly(
  data = umap,
  x = ~umap_1,
  y = ~umap_2,
  z = ~umap_3,
  color = ~seurat_clusters, # 按簇上色
  colors = cluster_colors,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
)

# 调整布局
plot_3d <- plot_3d %>% layout(
  scene = list(
    xaxis = list(title = "UMAP 1"),
    yaxis = list(title = "UMAP 2"),
    zaxis = list(title = "UMAP 3")
  )
)

# 保存为 PDF
library(webshot)
library(htmlwidgets)

# 保存为 HTML 再转为 PDF
html_file <- "all-Neuro-3D_scatterplot_colored.html"
pdf_file <- "3D_scatterplot_colored.pdf"

saveWidget(plot_3d, file = html_file)
webshot(html_file, pdf_file)


######拟时序空间#######

library(Seurat)
library(dplyr)
library(plotly)

# 1. 绘制降维图
DimPlot(af, pt.size = 0.2, group.by = "seurat_clusters", reduction = "umap")

# 2. 运行 3D UMAP 并提取嵌入
af <- RunUMAP(af, dims = 1:10, reduction = "pca", n.components = 3)
umap <- af@reductions$umap@cell.embeddings %>% as.data.frame()

# 3. 添加 Pseudotime 列
umap$Pseudotime <- af@meta.data$Pseudotime


# 提取 Monocle 的 Pseudotime 值
pseudotime <- as.data.frame(pData(cds)$Pseudotime)
rownames(pseudotime) <- rownames(pData(cds))
colnames(pseudotime) <- "Pseudotime"

seurat_object=sce

# 确保细胞名称一致
if (!all(rownames(pseudotime) %in% colnames(seurat_object))) {
  stop("细胞名称不匹配，请检查 Seurat 和 Monocle 对象是否来源于相同数据集")
}

# 将拟时间值添加到 Seurat 对象的元数据
seurat_object <- AddMetaData(seurat_object, metadata = pseudotime)


# 查看 Seurat 元数据中的 Pseudotime 列
head(seurat_object@meta.data)

# 可视化拟时间值分布
VlnPlot(seurat_object, features = "Pseudotime")
af=seurat_object



# 提取 UMAP 数据
umap <- af@reductions$umap@cell.embeddings %>% as.data.frame()

# 提取 Pseudotime 列并检查是否有空值
if ("Pseudotime" %in% colnames(af@meta.data)) {
  umap$Pseudotime <- af@meta.data$Pseudotime
} else {
  stop("af@meta.data does not contain a Pseudotime column.")
}

# 检查合并后列名
colnames(umap)



library(plotly)
library(htmlwidgets)
library(webshot)

# 绘制 3D 图，按 Pseudotime 上色
plot_3d <- plot_ly(
  data = umap,
  x = ~umap_1,
  y = ~umap_2,
  z = ~umap_3,
  color = ~Pseudotime,  # 使用 Pseudotime 映射颜色
  colors = colorRamp(c("#0000FF", "#7406d7", "#9d1da6", "#DC143C")),  # 定义颜色渐变
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
)

# 调整布局
plot_3d <- plot_3d %>% layout(
  scene = list(
    xaxis = list(title = "UMAP 1"),
    yaxis = list(title = "UMAP 2"),
    zaxis = list(title = "UMAP 3")
  )
)

# 保存为 HTML 文件
html_file <- "all-Neuro-3D_scatterplot_colored_by_pseudotime.html"
saveWidget(plot_3d, file = html_file)

# 安装 PhantomJS 后保存为 PDF
webshot::install_phantomjs()
pdf_file <- "3D_scatterplot_colored_by_pseudotime.pdf"
webshot(html_file, pdf_file)




########提取9，15############

setwd("D:\\王老师美国测序数据\\date")
load("25-01-06-6-cluster细胞注释.Rdata")


#####提取特定的new_cluster并合并#####
# 1. 创建新的 Seurat 对象，不包括 'Tcell', 'Ependymal cell', 'Micro/Macro' 细胞
# 检查 meta.data 中是否包含 new_cluster 列
head(scRNA_harmony@meta.data$new_cluster)
head(scRNA_harmony@meta.data)
# 将 new_cluster 列设为 Idents
Idents(scRNA_harmony) <- scRNA_harmony@meta.data$seurat_clusters


# 然后进行子集化，不包括 'Tcell', 'Ependymal cell', 'Micro/Macro' 细胞
scRNA_harmony_filtered <- subset(scRNA_harmony, 
                                 idents = c("9","15"))

colnames(scRNA_harmony@meta.data)

scRNA_harmony=scRNA_harmony_filtered

####聚类、umap/tsne降维降维####
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")

library(patchwork)

####调试dims####
for (i in c(5,10,15, 20, 25, 30)) {
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.2)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5) + labs(title = paste0("dims: ", i)))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)#assign() 函数将 plot_i 赋值给 plotname
  print(plot_i)
}
p = plot_5 + plot_10 + plot_15 + plot_20 + plot_25 + plot_30 +plot_layout( ncol = 3 )
ggsave( p, filename = "5-9dims.5-30.pdf", width = 15, height = 8 )
ggsave( plot_20, filename = "5-9dims.20.pdf", width = 5, height = 5 )


library(clustree)

set.seed(123)
# 再重新选择确定好的dim：10
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- FindClusters(
  object = scRNA_harmony,
  resolution = c( seq( 0.1, 1.5, 0.1) ) # 分辨率从 0.1,-1，间隔 0.1 一档
)
#删除一些RNA_snn_res
# scRNA_harmony@meta.data <- scRNA_harmony@meta.data[, -which(colnames(scRNA_harmony@meta.data) == "RNA_snn_res.1.5")]

# 可视化
clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.")
ggsave( filename = "5-9细分分群clustree.pdf", width = 12, height = 9 )


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
ggsave( p, filename = "5-9RNA_snn_res.0.1-1.5.pdf", width = 14, height = 15 )
p
p3 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5",
              label = T, repel = F, shuffle = T ,pt.size = 1,label.size = 5)
p3
ggsave(p3, filename = "5-9RNA_snn_res.0.5.pdf", width = 6, height = 6 )
# plot1 =DimPlot(scRNA_harmony, group.by = "RNA_snn_res.0.2",reduction = "umap",label = T,pt.size = 1,repel = F,label.size = 5)

#####正式降维参数选择####
#选择dim30，reso0.3
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.5)
table(scRNA_harmony@meta.data$seurat_clusters)


##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)

# 绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")#orig.ident
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = FALSE)#TRUE
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident")#celltype
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = FALSE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated
# 保存图片
ggsave("5-9_cells-umap_tsne_integrated1.pdf",umap_tsne_integrated,wi=10,he=10)

save(scRNA_harmony,file='9-15cluster—neure.Rdata')
DimPlot(scRNA_harmony)



markers = c("Snap25","Syp","Rbfox3","Snhg11",
            "Mbp","Mobp","Mog","Plp1",
            "Mpz","Pmp22","Prx",
            "Dcn","Col3a1","Igf2",
            "Aqp4","Atp1a2","Gja1","Slc1a2",
            "Flt1","Pecam1","Tek","Myl9","Pdgfrb",
            "Cspg4","Gpr17","Pdgfra",
            "Ctss","Itgam","Ptprc"
)

markers = c("Grp ","Tacr3 ","Reln",
            "Npy1r ","Nmur2 ","Car12",
            "Ntrk2 ","Pde11a ","Mybph",
            "Sstr2 ","Calb2 ","Bcl11b",
            "Gal ","Prkca ","Ankfn1",
            "Nos1 ","Pnoc ","Cdh8",
            "Prkcg ","Calb1 ","Ebf2",
            "Tac2 ","Necab1 ","Snca",
            "Nts ","Trpm3 ","Gbgt1",
            "Cck ","Ebb4 ","Ebf1",
            "Cbln2 ","Mgat5b ","Tshz2",
            "Adarb2 ","Rasa4 ","Dhrs3",
            "Maf ","Slc17a8 ","Parp8",
            "Cit ","Epha10 ","Kif5a",
            "Atrx ","Phf20l1 ","Kif21a",
            "Rsrp1  Stk39 ","Rab40c",
            "Gad2 ","Neto2 ","Htr2c",
            "Meis2  Sim1  Plaa",
            "Zfhx3  Ndst4  Grm8",
            "Zfhx4  Onecut  Vsx2",
            "Etl4 ","Fryl ","Abhd3",
            "Esrrg  Esrrb  Htr1f",
            "Tfap2b  Rab27b  Rarb",
            "Ssb ","Wdr47 ","Gls2",
            "Hapln1  Fhod3 ","Dgkg",
            "Lgr5 ","Foxp2 ","Slit2",
            "Tac1 ","Sox5 ","Otof",
            "Grpr ","Col5a2 ","Atp2b4",
            "Npy ","Ecel1 ","Qrfpr",
            "Man1a ","Iqgap2 ","Rreb1",
            "Mdga1 ","Megf11 ","Plcxd2",
            "Gabra1 ","Gabbr2 ","Cdh7",
            "Nfib ","Satb1 ","Rmst",
            "Lmx1b ","Pdyn ","Tcf4",
            "Col5a2 ","Qrfpr ","Penk",
            "Pitx2 ","Lyzl4 ","Chat",
            "Prph ","Isl1 ","Slit3",
            "Gbx1 ","Esr1 ","Mef2c",
            "Ptprz1 ","Zeb2 ","Zmat4",
            "Cdh3 ","Kcnip2 ","Sema3e",
            "Lhx1 ","Nxph1 ","Chrm3",
            "Nrgn ","Syt6 ","Fgf13",
            "Rorb ","Sorcs3 ","Col18a1")


markers = c("Olig3","Sim1",#V3
            "Shox2","Vsx2","Lhx3","Sox14","Lhx4",#V2a
            "Gata3","En1","Foxp2","Irx3","Dmrt3","Slc32a1",#V1/V2b
            "Foxp2","Foxp1","En1","Slc32a1",#V1
            "Dlk1","Slit2",#MN
            "Npy","Lhx1","Gad1","Gbx1","Kcnip2","Pax2","Rorb","Ptf1a",#dILA
            "Irx3","Dmrt3","Slc32a1",#dI6
            "Zic2","Pou4f1","Olig3","Lmx1b",#dI5/LB/V3
            "Cartpt","Prrxl1","Sst","Tlx3","Pou4f1","Tac1","Gal","Slc17a6","Cck","Mafa","Rora","Lbx1","Rorb","Lmx1b",#dI5/LB
            "Pou4f1","Bhlhe22","Lmx1b",#dI5
            "Tfap2b","Slc32a1","Slc6a5","Pdyn","Gal","Npy","Rorb","Mafb","Lhx1","Zic1","Mafa","Pax2","Zic2","Lbx1","Gad1","Ptf1a",#dI4
            "Isl1","Phox2a","Slc17a6",#dI3
            "Barhl2","Barhl1","Lhx2","Lhx9","Atoh1"#dI1
)
scRNA=scRNA_harmony
table(Idents(scRNA))

missing_genes <- markers[!markers %in% rownames(scRNA)]
print(missing_genes)
valid_markers <- markers[markers %in% rownames(scRNA)]
unique_markers <- unique(valid_markers)
length(valid_markers)
head(rownames(scRNA), 5)
DotPlot=DotPlot(scRNA, features = unique_markers) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#2166AC', '#67A9CF', '#D1E5F0', '#FDDBC7', '#EF8A62', '#B2182B'))
DotPlot
#neural_precursor_cells;
ggsave(filename="5-9-marker65.pdf",DotPlot, width = 150, height = 400, units = "mm")

####细胞比例#####


colour=c("#53b400","#f98b84","#a58aff","#00c094",
         "#fb61d7", "#c49a00","#00b6eb")



colour <- c("0" = "#faaca7",
            "10" = "#e1c066",
            "11" = "#ffa0e0",
            "13" = "#ddb0ff",
            "15" = "#66cbff",
            "17" = "#66d8db",
            "4" = "#66d8a0",
            "9" = "#b0ce66" )


#Idents(scRNA) <- scRNA$new_cluster

prop.table(table(Idents(scRNA)))

table(Idents(scRNA), scRNA$orig.ident)


Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("seurat_clusters", "sample", "ratio")


Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("seurat_clusters", "sample", "ratio")


ggplot(Cellratio, aes(x = sample, y = ratio, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "fill", colour = NA) +  # 去除线条
  scale_fill_manual(values = colour) +  # 使用自定义颜色
  labs(title = "Cell Type Proportions by Sample",
       x = "Sample",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


table(scRNA$orig.ident)
prop.table(table(Idents(scRNA)))  
table(Idents(scRNA), scRNA$orig.ident)  
Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)  
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("seurat_clusters","sample","ratio")



colourCount = length(unique(Cellratio$seurat_clusters)) 
ggplot(Cellratio) + 
  geom_bar(aes(x = sample, y = ratio, fill = seurat_clusters), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +  # 添加条形图层，根据细胞类型填充颜色
  theme_classic() + 
  labs(x = 'Sample', y = 'Ratio') +  
  coord_flip() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))+   # 设置面板边框属性
  scale_fill_manual(values = colour)  

p=ggplot(Cellratio) +
  geom_bar(aes(x =sample, y= ratio, fill = seurat_clusters),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = colour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p
# Save the plot as a PDF
ggsave("5-9cell_ratio_plot.pdf", plot = p, width = 4, height = 6)





library(Seurat)

setwd("D:\\王老师美国测序数据\\date")
load("9-15cluster—neure.Rdata")
#####12-18 神经元marker###### "Cdk4","Pcna","Eomes"#,"",
markers <- c("Gfap","Aldoc","Ebf1","Pou3f3",
             "Meis2","Rtn1","Map2","Erbb4",
             "Pde10a","Tshz2","Gria2","Rorb","Cux1","Thsd7a","Prox1","Sox1",
             "Notch1","Zfpm2","Cmip","Sox21","Foxp1","Foxp2","Plcb4","Rbfox3","Sox5",
             "Hopx","Egfr","Olig1","Olig2","Nkx2-2","Nkx6-2","Dbx2","Ascl1",
             "Sox2","Cdk4","Pcna","Mki67","Top2a","Sox2Gli3a","tdTomato")
#"Hoxa1","Hoxa2","Hoxa3","Hoxa4","Hoxa5","Hoxa6","Hoxa7","Hoxa9","Hoxa10","Hoxa11","Hoxa13","Hoxb1","Hoxb2","Hoxb3","Hoxb4","Hoxb5","Hoxb6","Hoxb7","Hoxb8","Hoxb9","Hoxb13","Hoxc4","Hoxc5","Hoxc6","Hoxc8","Hoxc9","Hoxc10","Hoxc11","Hoxc12","Hoxc13","Hoxd1","Hoxd3","Hoxd4","Hoxd8","Hoxd9","Hoxd10","Hoxd11","Hoxd12","Hoxd13")

# Check if scRNA is a Seurat object
scRNA_harmony=scRNA
# 检查当前的 cluster 标签
unique(scRNA_harmony$seurat_clusters)

# 如果有 NA 值，先移除它们
#scRNA_harmony <- subset(scRNA_harmony, subset = !is.na(seurat_clusters))

# 再次检查
unique(scRNA_harmony$seurat_clusters)

# 定义你想要的顺序
#desired_order <- c("0", "1",  "2", "3", "4")

# 检查所有指定的 levels 是否存在于数据中
missing_levels <- setdiff(desired_order, unique(scRNA_harmony$seurat_clusters))
if (length(missing_levels) > 0) {
  warning("以下 cluster 不存在于数据中: ", paste(missing_levels, collapse = ", "))
}

# 重新定义因子顺序
scRNA_harmony$seurat_clusters <- factor(scRNA_harmony$seurat_clusters, levels = desired_order)

# 验证顺序是否正确
levels(scRNA_harmony$seurat_clusters)

# 检查 markers 是否存在于数据中
valid_markers <- markers[markers %in% rownames(scRNA_harmony)]

# 生成 DotPlot
dotplot <- DotPlot(scRNA_harmony, features = valid_markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36')
  )
dotplot
# 显示图形
print(dotplot)

head(scRNA_harmony@meta.data$orig.ident)
table(scRNA_harmony@meta.data$orig.ident)
# 合并 cluster 和 orig.ident 信息
scRNA_harmony$cluster_group <- paste0(scRNA_harmony$seurat_clusters, "_", scRNA_harmony$orig.ident)

# 重新设置顺序（保证 cluster 顺序不乱）
scRNA_harmony$cluster_group <- factor(
  scRNA_harmony$cluster_group,
  levels = unlist(lapply(desired_order, function(cl) paste0(cl, "_", c("sgrna", "sorna"))))
)

# 生成 DotPlot
dotplot <- DotPlot(scRNA_harmony, features = valid_markers, group.by = "cluster_group") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36')
  )
dotplot


ggsave(filename = "8-28-5-9-GLI3-all-cluster-marker.pdf", plot = dotplot, width = 115, height = 210, units = "mm")
#####差异分析#########

# 差异分析
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

#2.对计算好的每cluster的marker基因进行筛选
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
View(top10)
write.csv(top10,"5cluster_top10.csv",row.names = T)
write.csv(top20,"5cluster_top20.csv",row.names = T)
write.csv(top50,"5cluster_top50.csv",row.names = T)

top10 <- as.data.frame(top10)  # 将markers对象转换为数据框
markerdata <- ScaleData(scRNA_harmony, features = as.character(unique(top10$gene)), assay = "RNA")  # 对特定基因进行标准化
# 常规热图
DoHeatmap(markerdata,  # 使用DoHeatmap函数绘制热图
          features = as.character(unique(top10$gene)),  # 指定要显示的特征
          group.by = "seurat_clusters",  # 指定按照免疫细胞类型分组celltype
          assay = 'RNA')  # 指定分析的分析类型为RNA


#修改颜色
DoHeatmap(markerdata,  # 再次绘制热图，附加自定义颜色和渐变色
          features = as.character(unique(top10$gene)),
          group.by = "seurat_clusters",#celltype
          assay = 'RNA',
          group.colors = c("#B288AB","#D56B19","#85AD9B","#747D90","#D93B45" ,"#8DBCCA",
                           "#504646","#C89B60","#FAEACE","#F8B9AB" ,"#F9F4F5")) +
  scale_fill_gradientn(colors = c("#008D83","grey","firebrick3"))


ggsave(filename = "5-9DoHeatmap热图1.pdf", width = 12, height = 12 )



######ClusterGVis绘制单细胞亚群marker基因###########
library(Seurat)
library(dplyr)
library(ggplot2)
library(SCP)
library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(org.Mm.eg.db)#org.Mm.eg.db
sce=scRNA_harmony
table(sce$seurat_clusters)#celltype
Idents(sce) <- 'seurat_clusters'#celltype

DimPlot(sce,label = T)


markers<-FindAllMarkers(sce, only.pos=TRUE,
                        min.pct=0.25,
                        logfc.threshold=0.25)#get top 20genes
markers<-deg
markers<-markers%>%
  dplyr::group_by(cluster)%>%
  dplyr::top_n(n=20,wt=avg_log2FC)


head(markers)

data<-prepareDataFromscRNA(object=sce,
                           diffData=markers,
                           showAverage=TRUE)
str(data)

enrich<-enrichCluster(object=data,
                      OrgDb=org.Mm.eg.db,
                      type="BP",
                      organism="mmu",
                      pvalueCutoff=0.05,
                      topn=5,
                      seed=123)

#check
head(enrich)

#add gene name
markGenes=unique(markers$gene)[sample(1:length(unique(markers$gene)),35,
                                      replace=F)]

#line plot
visCluster(object=data,
           plot.type="line")


visCluster(
  object = data,
  plot.type  = "heatmap",
  column_names_rot = 35,
  markGenes = markGenes,
  cluster.order  = c(1:7)
)

pdf("5-9热图带富集4.pdf",width = 16,height = 18)
visCluster(
  object          = data,
  plot.type        = "both",
  column_names_rot = 25,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side   = "left",
  annoTerm.data    = enrich,
  line.side        = "left",
  cluster.order    = c(1:25),
  # 假设enrich有41个条目
  go.col = rep(jjAnno::useMyCol("stallion", n = 5), each = 5)[1:25],
  add.bar = TRUE
)
dev.off()

#####SCP差异基因分析######
library(SCP)
library(BiocParallel)
register(MulticoreParam(workers = 12, progressbar = TRUE))
data("sce")
print(sce)
setwd("D:\\王老师美国测序数据\\date")
load("9-15cluster—neure.Rdata")

sce <- RunDEtest(srt = sce, group_by = "seurat_clusters", fc.threshold = 1, only.pos = FALSE)#celltype
VolcanoPlot(srt = sce, group_by = "seurat_clusters")


DEGs <- sce@tools$DEtest_seurat_clusters$AllMarkers_wilcox


DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
library(Seurat)
library(SingleCellExperiment)

names(sce@assays)
sce <- AnnotateFeatures(sce, species = "Mus_musculus", db = c("TF", "CSPA"))#Mus_musculus,Homo_sapiens
ht <- FeatureHeatmap(
  srt = sce, group.by = "seurat_clusters", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Mus_musculus", db = c("GO_BP", "KEGG"), anno_terms = TRUE,#Region
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)
print(ht$plot)
ggsave("5-9-1heatmap1.pdf", plot = ht$plot, width = 24, height = 12, units = "in")

library(patchwork)
library(ggplot2)
library(ggplot2)
library(patchwork)

# 定义新的基因列表，选择具有更高特异性的标记基因
g <- c("Erbb4",
       "Pde10a","Tshz2","Gria2","Rorb","Cux1","Thsd7a","Prox1","Sox1",
       "Notch1","Zfpm2","Cmip","Sox21","Foxp1","Foxp2","Plcb4"
)
#化学门控通道基因

# 生成一个存储多个 ggplot 对象的列表
p_merge <- list()

# 循环绘制每个基因的图
for (i in 1:length(g)) {
  # 打印当前基因
  print(g[i])
  
  # 画图
  p_merge[[i]] <- FeaturePlot(sce, features = g[i], order = TRUE, pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#e0f2f7", "#b3d8e6", "#80bed9", "#4aa3c8", "#1f8abf", "#1570a6", "#0f5688"),
      limits = c(0, 4),
      name = "Expression"# 所有图表使用相同的图例名称
    ) +
    theme_void() +
    labs(title = g[i]) +  # 添加基因名作为标题
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"# 每个单独图表不显示图例
    )
}

# 创建一个单独的图例
legend_plot <- FeaturePlot(sce, features = g[1], order = TRUE) +
  scale_color_gradientn(
    colors = c("#e0f2f7", "#b3d8e6", "#80bed9", "#4aa3c8", "#1f8abf", "#1570a6", "#0f5688"),
    limits = c(0, 4),
    name = "Expression"
  ) +
  theme(legend.position = "right")


# 提取图例
legend <- cowplot::get_legend(legend_plot)

# 使用 patchwork 的 wrap_plots() 函数组合图表
plot_grid <- wrap_plots(p_merge, ncol = 3)  # 3x3 布局

# 添加图例到网格
combined_plot <- cowplot::plot_grid(
  plot_grid,
  legend,
  rel_widths = c(1, 0.2)
)

# 显示合并后的图
combined_plot

# 或者使用基础图形设备
pdf("8-28-all_plot.pdf", width=12, height=9)
print(combined_plot)
dev.off()


save(sce,file='9-15-4cluster—neure.Rdata')

