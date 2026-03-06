rm(list=ls())
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

# load 10x data 
C6_data <- Read10X(data.dir = "/mnt/data5/disk/wupeng/A549_lineage/2.scRNA_pepline/1.cellranger_result/C6_count/outs/filtered_feature_bc_matrix/")
E8_data <- Read10X(data.dir = "/mnt/data5/disk/wupeng/A549_lineage/2.scRNA_pepline/1.cellranger_result/E8_count/outs/filtered_feature_bc_matrix/")
G5_data <- Read10X(data.dir = "/mnt/data5/disk/wupeng/A549_lineage/2.scRNA_pepline/1.cellranger_result/G5_count/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data)
sc_C6 <- CreateSeuratObject(C6_data, project = "A549_C6", min.cells = 3, min.features = 200)
sc_E8 <- CreateSeuratObject(E8_data, project = "A549_E8", min.cells = 3, min.features = 200)
sc_G5 <- CreateSeuratObject(G5_data, project = "A549_G5", min.cells = 3, min.features = 200)

sc_C6[['percent.mt']] <- PercentageFeatureSet(sc_C6, pattern = 'MT-')
sc_E8[['percent.mt']] <- PercentageFeatureSet(sc_E8, pattern = 'MT-')
sc_G5[['percent.mt']] <- PercentageFeatureSet(sc_G5, pattern = 'MT-')


VlnPlot(object = sc_C6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = sc_E8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = sc_G5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sc_C6, feature1 = "nCount_RNA", feature2 = "percent.mt")+ NoLegend()+labs(x="nUMI",y="percent.mt") #-0.20
plot2 <- FeatureScatter(sc_C6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +labs(x="nUMI", y="nGene") #0.92
plot1+plot2
plot1 <- FeatureScatter(sc_E8, feature1 = "nCount_RNA", feature2 = "percent.mt")+ NoLegend()+labs(x="nUMI",y="percent.mt") #-0.22
plot2 <- FeatureScatter(sc_E8, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +labs(x="nUMI", y="nGene") #0.91
plot1+plot2
plot1 <- FeatureScatter(sc_G5, feature1 = "nCount_RNA", feature2 = "percent.mt")+ NoLegend()+labs(x="nUMI",y="percent.mt") #-0.24
plot2 <- FeatureScatter(sc_G5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +labs(x="nUMI", y="nGene") #0.94
plot1+plot2

##filter MT genes
C6.filterMT <- subset(sc_C6, subset = nFeature_RNA > 200 & percent.mt < 25) #9440
C6.filterMT <- RenameCells(C6.filterMT,add.cell.id = "A549_C6")
dim(C6.filterMT)
E8.filterMT <- subset(sc_E8, subset = nFeature_RNA > 200 & percent.mt < 25) #3533
E8.filterMT <- RenameCells(E8.filterMT,add.cell.id = "A549_E8")
dim(E8.filterMT)
G5.filterMT <- subset(sc_G5, subset = nFeature_RNA > 200 & percent.mt < 25) #5400
G5.filterMT <- RenameCells(G5.filterMT,add.cell.id = "A549_G5")
dim(G5.filterMT)

##load Seurat V5
C6.V5 <- readRDS("/mnt/data5/disk/wupeng/A549_lineage/2.scRNA_pepline/2.filter_mt_doublet/filter_0.25mt_doublet/A549_C6_seuratObj.Rds")
dim(C6.V5@meta.data) #8727
E8.V5 <- readRDS("/mnt/data5/disk/wupeng/A549_lineage/2.scRNA_pepline/2.filter_mt_doublet/filter_0.25mt_doublet/A549_E8_seuratObj.Rds")
dim(E8.V5@meta.data) #3433
G5.V5 <- readRDS("/mnt/data5/disk/wupeng/A549_lineage/2.scRNA_pepline/2.filter_mt_doublet/filter_0.25mt_doublet/A549_G5_seuratObj.Rds")
dim(G5.V5@meta.data) #5167

cells_in_C6.V5 <- rownames(C6.V5@meta.data)
C6.filterMT.V4 <- subset(C6.filterMT, cells = cells_in_C6.V5)
dim(C6.filterMT.V4)

cells_in_E8.V5 <- rownames(E8.V5@meta.data)
E8.filterMT.V4 <- subset(E8.filterMT, cells = cells_in_E8.V5)
dim(E8.filterMT.V4)

cells_in_G5.V5 <- rownames(G5.V5@meta.data)
G5.filterMT.V4 <- subset(G5.filterMT, cells = cells_in_G5.V5)
dim(G5.filterMT.V4)


Seurat.list <- list(C6.filterMT.V4, E8.filterMT.V4, G5.filterMT.V4)
all.Seuratobj <- Reduce(function(x,y) merge(x,y), Seurat.list)
saveRDS(all.Seuratobj, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/all.rawSeuratobj.V4.Rds")




rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(SingleCellExperiment)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

#============Run===============================================================================
#==1 subset same number cells for HSC and hESC samples
all.Seuratobj <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/all.rawSeuratobj.V4.Rds")
table(all.Seuratobj$orig.ident)
scobj <- NormalizeData(all.Seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
#==2 Find highly variable features
scobj <- FindVariableFeatures(scobj, selection.method = "vst", verbose = T, nfeatures = 3000)
scobj <- CellCycleScoring(scobj, s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)
#==3 ScaleData
scobj <- ScaleData(scobj, features = rownames(scobj))
#scobj <- RunPCA(scobj, features = c(cc.genes$s.genes,cc.genes$g2m.genes))
#DimPlot(scobj,reduction = "pca")
#DimPlot(scobj,reduction = "pca",split.by='Phase',ncol=1)
#DimPlot(scobj,reduction = "pca",group.by="orig.ident")
#DimPlot(scobj,reduction = "pca",split.by='orig.ident',ncol=2)

#==4 RUN PCA
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),verbose = F, npcs = 100)
ElbowPlot(scobj)
ElbowPlot(scobj, reduction = "pca",ndims=15)
DimPlot(scobj, reduction = "pca",group.by="orig.ident")

scobj <- RunUMAP(scobj, reduction = "pca",dims=1:15)
DimPlot(scobj, reduction = "umap",group.by="orig.ident")

#==5 harmony
scobj.batch <- RunHarmony(scobj, group.by.vars="orig.ident")
#scobj.batch <- RunHarmony(scobj.ref,group.by.vars="orig.ident", lambda=NULL)
scobj.batch <- RunUMAP(scobj.batch, reduction = "harmony",dims=1:15) #, reduction.name = "umap")
DimPlot(scobj.batch, reduction="umap",group.by="orig.ident",pt.size=0.5)

#==6 cluster
scobj.batch <- FindNeighbors(scobj.batch, reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.2) #10
DimPlot(scobj.batch, reduction="umap", label = TRUE,pt.size=0.5)
saveRDS(scobj.batch, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/scobj.V4.Rds")

#==7 find marker
rm(list=ls())
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/scobj.V4.Rds")
sce$celltype <- plyr::mapvalues(x=sce$seurat_clusters,
                                from = c("2","3","0","5","4","1"),
                                to = c(paste0("C",seq(0,5))))
sce$celltype <- factor(sce$celltype, levels = c(paste0("C",seq(0,5))))
Idents(sce) <- "celltype"

all.marker.1 <- FindAllMarkers(sce, only.pos=TRUE,min.pct=0.1,logfc.threshold = 0.25,test.use = "wilcox")
top10 <- all.marker.1 %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% ungroup()
top50 <- all.marker.1 %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 50) %>% ungroup()
top100 <- all.marker.1 %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 100) %>% ungroup()
#save(all.marker.1, top10, top50, top100,
#     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/Marker.geneList.Rda")
write.csv(top10,file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/Marker.top10.csv")


top10_log2FC_0.5 <- all.marker.1 %>% group_by(cluster) %>% dplyr::filter(avg_log2FC >= 0.5) %>% 
  arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% ungroup()
#top10 %>% filter(cluster==5) %>% select(gene) %>% pull()
top10_log2FC_1.0 <- all.marker.1 %>% group_by(cluster) %>% dplyr::filter(avg_log2FC >= 1.0) %>% 
  arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% ungroup()

clusterCols <- c("#B2DF8A","#FE9999","#A4CEDE","#1F78B4","#3DA42D","#BBBFE0")
#display.brewer.pal(11,"RdBu")
#brewer.pal(11,"RdBu")

p1 <- DoHeatmap(sce, features = top10$gene, group.by="celltype",group.colors=clusterCols,draw.lines=T) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n=11,name="RdBu")),na.value="white") +
  #scale_fill_gradientn(colors = colorRampPalette(c("#92C5DE", "#F4A582"))(10)) +
  theme(axis.text.x=element_blank(),axis.text.y=element_text(size=5,color="black"),
        axis.line = element_blank(),axis.ticks = element_blank())
p2 <- DoHeatmap(sce, features = top10_log2FC_0.5$gene, group.by="celltype",group.colors=clusterCols,draw.lines=T) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n=11,name="RdBu")),na.value="white") +
  #scale_fill_gradientn(colors = colorRampPalette(c("#92C5DE", "#F4A582"))(10)) +
  theme(axis.text.x=element_blank(),axis.text.y=element_text(size=5,color="black"),
        axis.line = element_blank(),axis.ticks = element_blank())
#p1|p2

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/0.2.top10_MarkerGenes.Heatmap.pdf",width=5,height=12)
pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/Fig.4-28.pdf",width=6,height=6)
p1
dev.off()

##tree info
rm(list=ls())
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/scobj.V4.Rds")
sce$celltype <- plyr::mapvalues(x=sce$seurat_clusters,
                                from = c("2","3","0","5","4","1"),
                                to = c(paste0("C",seq(0,5))))
###检测每种细胞类型的数目,与tree信息对应
cell <- sce@meta.data %>% dplyr::select(c(orig.ident, seurat_clusters, celltype)) %>%
  mutate(BC=rownames(.)) %>% mutate(BC=gsub("-1","",.$BC)) %>%
  mutate(BC=gsub("A549_C6|A549_E8|A549_G5", "BC", .$BC)) %>% as.data.frame()
unique(cell$orig.ident)
cell %>% select(orig.ident,celltype) %>% group_by(orig.ident) %>% table()

allele_path <- "/mnt/data5/disk/wupeng/A549_lineage/2024_5_30_KCA/filter_100indel/seuratOBJ_cell_type_info/"
C6_treeInfo <- readRDS(paste0(allele_path, "A549_C6","_all_cell_info.RDS"))
cell_C6 <- cell %>% filter(orig.ident=="A549_C6")
C6_treeInfo <- C6_treeInfo %>% mutate(celltype=cell_C6$celltype[match(C6_treeInfo$BC,cell_C6$BC)],
                                      seurat_clusters=cell_C6$seurat_clusters[match(C6_treeInfo$BC,cell_C6$BC)])
table(C6_treeInfo$celltype)
C6_treeInfo_single <- C6_treeInfo %>% filter(cell_num==1)
table(C6_treeInfo_single$celltype)
#
E8_treeInfo <- readRDS(paste0(allele_path, "A549_E8","_all_cell_info.RDS"))
cell_E8 <- cell %>% filter(orig.ident=="A549_E8")
E8_treeInfo <- E8_treeInfo %>% mutate(celltype=cell_E8$celltype[match(E8_treeInfo$BC,cell_E8$BC)],
                                      seurat_clusters=cell_E8$seurat_clusters[match(E8_treeInfo$BC,cell_E8$BC)])
table(E8_treeInfo$celltype)
E8_treeInfo_single <- E8_treeInfo %>% filter(cell_num==1)
table(E8_treeInfo_single$celltype)
#
G5_treeInfo <- readRDS(paste0(allele_path, "A549_G5","_all_cell_info.RDS"))
cell_G5 <- cell %>% filter(orig.ident=="A549_G5")
G5_treeInfo <- G5_treeInfo %>% mutate(celltype=cell_G5$celltype[match(G5_treeInfo$BC,cell_G5$BC)],
                                      seurat_clusters=cell_G5$seurat_clusters[match(G5_treeInfo$BC,cell_G5$BC)])
table(G5_treeInfo$celltype)
G5_treeInfo_single <- G5_treeInfo %>% filter(cell_num==1)
table(G5_treeInfo_single$celltype)

save(C6_treeInfo, E8_treeInfo, G5_treeInfo, 
     file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/treeInfo.Rda")




#==plot===========================================================================================
rm(list=ls())
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/scobj.V4.Rds")
sce$celltype <- plyr::mapvalues(x=sce$seurat_clusters,
                                from = c("2","3","0","5","4","1"),
                                to = c(paste0("C",seq(0,5))))
sce$celltype <- factor(sce$celltype, levels = c(paste0("C",seq(0,5))))
#1 UMAP sample
sample <- c("A549_C6","A549_E8","A549_G5")
color <- c("#1F78B4","#A6CEE3","#B2DF8A")
#color <- c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#FB9A99","#6A3D9A","#FDBF6F")
color <- setNames(color, sample)

sce$orig.ident <- factor(sce$orig.ident, levels = sample)
p1 <- DimPlot(sce,reduction="umap",group.by="orig.ident",pt.size=0.1,cols = color) +
  labs(x="UMAP1",y="UMAP2",title="")

#2 UMAP cluster
#clusterCols <- c("#B2DF8A","#3DA42D","#FE9999","#A4CEDE","#FDBE6F","#FF9116","#BBBFE0","#673A95","#B25A2C")
clusterCols <- c("#B2DF8A","#FE9999","#A4CEDE","#1F78B4","#3DA42D","#BBBFE0")
p2 <- DimPlot(sce, reduction="umap",label = T, cols=clusterCols, pt.size=0.1, group.by="celltype")+
  labs(x="UMAP1",y="UMAP2",title="")

#3  Cluster distribution per sample
cluster <- sce@meta.data[,c(1,11)]
res1 <- cluster %>% group_by(orig.ident) %>% dplyr::count(celltype) %>% ungroup()
res2 <- res1 %>% dplyr::group_by(orig.ident) %>% dplyr::summarise(n_sum = sum(n)) %>% ungroup()
res <- merge(res1,res2,by='orig.ident') %>% mutate(ratio=n/n_sum) %>% select(1,2,5)
res$orig.ident <- factor(res$orig.ident, levels=c("A549_C6","A549_E8","A549_G5"))

p3 <- ggplot(res,aes(x=orig.ident,y=ratio,fill=celltype))+geom_bar(stat="identity", width=0.6, position="fill",alpha=0.8)+
  labs(x="",y="")+scale_fill_manual(values=clusterCols)+ #theme_classic()+ 
  #scale_y_continuous(expand=c(0,0))+
  theme(panel.border=element_rect(fill='transparent',color='black'),panel.background=element_blank(),
        axis.ticks=element_line(color="black"),axis.ticks.length=unit(0.13,"cm"),
        axis.title.y = element_text(size=12,color="black"),axis.text.x=element_text(size=10,color='black',angle=30,hjust=1,vjust=1),
        axis.text.y=element_text(size=10,color="black"))


pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/0.2.umap.pdf",width=11.5,height=5)
p1|p2
dev.off()

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/0.2.distribution.pdf",width=3.5,height=5)
p3
dev.off()








