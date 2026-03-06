rm(list=ls())
## ==================== import some librarys =====================================================
suppressMessages({
  library(SOT)
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

## ============Run================================================================================
#====== 0 load and merge Seurat V4 object =====================
B7 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scHSC_B7_seuratObj.Rds")
F12 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scHSC_F12_seuratObj.Rds")
E12 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scHSC_E12_seuratObj.Rds")
F8 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scHSC_F8_seuratObj.Rds")
A7 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scH9_A7_seuratObj.Rds")
C9 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scH9_C9_seuratObj.Rds")
C12 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scH9_C12_seuratObj.Rds")
D3 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/seuratObj_Mt0.25/scH9_D3_seuratObj.Rds")

Seurat.list <- list(B7, E12, F12, A7, C9, C12, D3)
all.Seuratobj<-Reduce(function(x,y) merge(x,y), Seurat.list)
saveRDS(all.Seuratobj, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/all.rawSeuratobj.V4.Rds")

#====== 1 Mapping and annotating query datasets ===============
## Step 1 subset same number cells for HSC and hESC samples
all.Seuratobj <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/all.rawSeuratobj.V4.Rds")

allcell <- all.Seuratobj@meta.data %>% mutate(cell=rownames(.))
cell.ref.1 <- allcell %>% filter(orig.ident %in% c("10X_HSC_B7", "10X_HSC_E12", "10X_HSC_F12"))
cell.ref.2 <- allcell %>% filter(!orig.ident %in% c("10X_HSC_B7", "10X_HSC_E12", "10X_HSC_F12")) %>%
  sample_n(size = nrow(cell.ref.1), replace = FALSE)
table(cell.ref.2$orig.ident)
cell.ref <- rbind(cell.ref.1, cell.ref.2)
table(cell.ref$orig.ident)
##
cell.query <- allcell %>% filter(!cell %in% cell.ref$cell)
table(cell.query$orig.ident)

## Step 2 Cell type classification using an integrated reference
scobj.ref <- subset(all.Seuratobj, cells=cell.ref$cell)
table(scobj.ref@meta.data$orig.ident)
#save(scobj.ref, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/ref.Seuratobj.V4.Rda")

#==2.1 get hvg.genes list
#The data manipulation of SOT
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/ref.Seuratobj.V4.Rda")
scobj.ref <- NormalizeData(scobj.ref, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- as.SingleCellExperiment(scobj.ref)
range(sce@assays@data$logcounts) #0.00000 8.36534
#Filter low expression gene
t <- as.data.frame(sce@assays@data$logcounts["ECM1",])
length(which(t[,1]!=0))/length(t[,1])
max(sce@assays@data$logcounts["PDGFRB",]) #2.588912
max(sce@assays@data$logcounts["COL1A1",]) #4.071364
max(sce@assays@data$logcounts["ACTA2",])  #3.93428
max(sce@assays@data$logcounts["ALCAM",])  #3.452557
max(sce@assays@data$logcounts["PCDH7",])  #3.363185
max(sce@assays@data$logcounts["CCL5",])   #6.053652
max(sce@assays@data$logcounts["COL3A1",]) #4.081807

sce <- FilterLow(sce, minexp = 3, mincells = 3, datatype = "logcounts") ## minexp factor is key
#Find high variable genes
sce <- FindHVGs(sce, datatype = "logcounts", thr.bio = 0, thr.FDR = 0.1)
t2 <- as.data.frame(rowData(sce)[rowData(sce)$genes.use,]) #914

hvg.genes <- rowData(sce)$symbol[rowData(sce)$genes.use] %>% as.character()
hsc_marker <- c("DES","ALCAM","ACTA2","COL1A1","LRAT","RELN","PCDH7","PDGFRB","ECM1","POU5F1","CCL5","COL3A1")
hsc_marker %in% hvg.genes
hvg.genes.add <- c(hvg.genes,"DES","LRAT","RELN","PDGFRB","ECM1")
#writeLines(hvg.genes.add,"/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.hvg.genes.add.txt")

#==2.3 ScaleData
# Find highly variable features
scobj.ref <- FindVariableFeatures(scobj.ref, selection.method = "vst", verbose = T, nfeatures = 3000)
# cell cycle scoring
scobj.ref <- CellCycleScoring(scobj.ref, s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)
# scale data
scobj.ref <- ScaleData(scobj.ref, features = rownames(scobj.ref))

#==2.4 RUN PCA for seurat object by hvg.genes of SOT
scobj.ref <- RunPCA(scobj.ref, features=hvg.genes.add, verbose = F, npcs = 100)
ElbowPlot(scobj.ref)
ElbowPlot(scobj.ref, reduction = "pca",ndims=15)
DimPlot(scobj.ref, reduction = "pca",group.by="orig.ident")

scobj.ref <- RunUMAP(scobj.ref, reduction = "pca",dims=1:15)
DimPlot(scobj.ref, reduction = "umap",group.by="orig.ident")

#==2.5 harmony
scobj.batch <- RunHarmony(scobj.ref, group.by.vars="orig.ident")
#scobj.batch <- RunHarmony(scobj.ref,group.by.vars="orig.ident", lambda=NULL)
scobj.batch <- RunUMAP(scobj.batch, reduction = "harmony",dims=1:15) #, reduction.name = "umap")
DimPlot(scobj.batch, reduction="umap",group.by="orig.ident",pt.size=0.5)

#==2.6 cluster
scobj.batch <- FindNeighbors(scobj.batch, reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.2) #10
DimPlot(scobj.batch, reduction="umap", label = TRUE,pt.size=0.5)

###number of each cell type
cell <- scobj.batch@meta.data %>% dplyr::select(c(orig.ident,seurat_clusters))
cell$BC <- rownames(cell)
cell$BC <- gsub("-1","",cell$BC)
cell$BC <- gsub(".*SC_", "BC_", cell$BC)
cell$orig.ident <- gsub("10X_","",cell$orig.ident)
cell$orig.ident <- gsub("H9_","hESC_",cell$orig.ident)
unique(cell$orig.ident)

sample.name="hESC_D3"
sample.name="HSC_B7"
sample.name="HSC_E12"
sample.name="HSC_F12" #"hESC_D3"
allelesInfo_path = paste0("/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/",
                          sample.name, "_B/", sample.name, ".AllelesInfo.csv")
alleinfo <- read.csv(allelesInfo_path)
celltype <- subset(cell,cell$orig.ident==sample.name)
#celltype$celltype <- as.character(celltype$celltype)
alleinfo.celltype <- merge(alleinfo,celltype,by="BC") %>% dplyr::select(nodeLabel,BC,seurat_clusters,cellNum)
One.celltype <- alleinfo.celltype%>%filter(cellNum==1) 
table(One.celltype$seurat_clusters)
table(alleinfo.celltype$seurat_clusters)
table(celltype$seurat_clusters)

#==2.7 define celltype base on other marker gene list ==================================================
D4_gene <- read.csv("/mnt/data5/disk/yangwj/3.3.SOT/DEGs/HSC_4d.CSV") 
D6_gene <- read.csv("/mnt/data5/disk/yangwj/3.3.SOT/DEGs/HSC_6d.CSV")
D8_gene <- read.csv("/mnt/data5/disk/yangwj/3.3.SOT/DEGs/HSC_8d.CSV")
D12_gene <- read.csv("/mnt/data5/disk/yangwj/3.3.SOT/DEGs/HSC_12d.CSV")
D0.array.gene.2 <- read.csv("/mnt/data5/disk/yangwj/3.3.SOT/DEGs/day0.array.undifferentiated.iPSC.padj0.05.csv",row.names = 1)
D12.array.gene.2 <- read.csv("/mnt/data5/disk/yangwj/3.3.SOT/DEGs/day12.array.iPSC-HSC.padj0.05.csv",row.names = 1)

sce <- scobj.batch
Idents(sce) <- "seurat_clusters"
#== get Score========
#==D4
genelist <- list(D4_gene$D5_names)
sce <- Seurat::AddModuleScore(object=sce, features=genelist,ctrl = 100,name="D4.sc")
colnames(sce@meta.data)[12] <- "D4_sc_Score"
D4_sc_Score <- FetchData(sce,vars = c("seurat_clusters","D4_sc_Score"))
#==D6
genelist <- list(D6_gene$D7_names)
sce <- Seurat::AddModuleScore(object=sce, features=genelist,ctrl = 100,name="D6.sc")
colnames(sce@meta.data)[13] <- "D6_sc_Score"
D6_sc_Score <- FetchData(sce,vars = c("seurat_clusters","D6_sc_Score"))
#==D8
genelist <- list(D8_gene$D9_names)
sce <- Seurat::AddModuleScore(object=sce, features=genelist,ctrl = 100,name="D8.sc")
colnames(sce@meta.data)[14] <- "D8_sc_Score"
D8_sc_Score <- FetchData(sce,vars = c("seurat_clusters","D8_sc_Score"))
#==D12
genelist <- list(D12_gene$D13_names)
sce <- Seurat::AddModuleScore(object=sce, features=genelist,ctrl = 100,name="D12.sc")
colnames(sce@meta.data)[15] <- "D12_sc_Score"
D12_sc_Score <- FetchData(sce,vars = c("seurat_clusters","D12_sc_Score"))
#==D0_arrray
genelist <- list(D0.array.gene.2$gene)
sce <- Seurat::AddModuleScore(object=sce, features=genelist,ctrl = 100,name="D0.array")
colnames(sce@meta.data)[16] <- "D0_array_Score"
D0_array_Score <- FetchData(sce,vars = c("seurat_clusters","D0_array_Score"))
#==D12_arrray
genelist <- list(D12.array.gene.2$gene)
sce <- Seurat::AddModuleScore(object=sce, features=genelist,ctrl = 100,name="D12.array")
colnames(sce@meta.data)[17] <- "D12_array_Score"
D12_array_Score <- FetchData(sce,vars = c("seurat_clusters","D12_array_Score"))

p1 <- DotPlot(sce,features=c("D0_array_Score","D4_sc_Score","D6_sc_Score","D8_sc_Score","D12_sc_Score","D12_array_Score"),
              cols = "RdYlBu",dot.scale=7)+labs(x="",y="")+ theme(axis.text.x=element_text(angle=30,hjust=1))
p2 <- DotPlot(sce,features=c("D0_array_Score","D4_sc_Score","D8_sc_Score","D12_sc_Score","D12_array_Score"),
              cols = "RdYlBu",dot.scale=7)+labs(x="",y="")+ theme(axis.text.x=element_text(angle=30,hjust=1))
p1+p2

##rename clusters
Idents(sce, cells = WhichCells(sce, idents = c(7, 9))) <- 7
sce@meta.data$population <- Idents(sce)
DimPlot(sce, reduction="umap", label = TRUE,pt.size=0.5, group.by = "population")
sce$celltype <- plyr::mapvalues(x=sce$population,
                                from = c("0","6","5","3","7","4","1","2","8"),
                                to = c(paste0("C",seq(0,8))))
sce$celltype <- factor(sce$celltype,levels=c(paste0("C",seq(0,8))))
saveRDS(sce, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")


## Step3 RUN PCA on the remaining hESC cells using hvg.genes identified by SOT. Then, ref and query were merged and re-grouped.
allcell <- all.Seuratobj@meta.data %>% mutate(cell=rownames(.))
cell.query <- allcell %>% filter(!cell %in% rownames(scobj.batch@meta.data))

scobj.query <- subset(all.Seuratobj, cells=cell.query$cell)
table(scobj.query@meta.data$orig.ident)

scobj.query <- NormalizeData(scobj.query, normalization.method = "LogNormalize", scale.factor = 10000)
scobj.query <- FindVariableFeatures(scobj.query, selection.method = "vst", verbose = T, nfeatures = 3000)
scobj.query <- CellCycleScoring(scobj.query, s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)
scobj.query <- ScaleData(scobj.query, features = rownames(scobj.query))

scobj.query <- RunPCA(scobj.query, features=hvg.genes.add, verbose = F, npcs = 100)
ElbowPlot(scobj.query)
ElbowPlot(scobj.query, reduction = "pca",ndims=15)
DimPlot(scobj.query, reduction = "pca",group.by="orig.ident")

scobj.query<- RunUMAP(scobj.query, reduction = "pca",dims=1:15)
DimPlot(scobj.query, reduction = "umap",group.by="orig.ident")

scobj.query.batch <- RunHarmony(scobj.query, group.by.vars="orig.ident")
scobj.query.batch <- RunUMAP(scobj.query.batch, reduction = "harmony",dims=1:15)
DimPlot(scobj.query.batch, reduction="umap",group.by="orig.ident",pt.size=0.5)
#cluster
scobj.query.batch <- FindNeighbors(scobj.query.batch, reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.2) #5
DimPlot(scobj.query.batch, reduction="umap", label = TRUE,pt.size=0.5)
##rename clusters
Idents(scobj.query.batch, cells = WhichCells(scobj.query.batch, idents = c(0, 1))) <- 0
Idents(scobj.query.batch, cells = WhichCells(scobj.query.batch, idents = 3)) <- 5
Idents(scobj.query.batch, cells = WhichCells(scobj.query.batch, idents = 2)) <- 3
Idents(scobj.query.batch, cells = WhichCells(scobj.query.batch, idents = 4)) <- 6
scobj.query.batch@meta.data$population <- Idents(scobj.query.batch)
DimPlot(scobj.query.batch, reduction="umap", label = TRUE,pt.size=0.5, group.by = "population")
scobj.query.batch$celltype <- plyr::mapvalues(x=scobj.query.batch$population,
                                              from = c("0","6","5","3"),
                                              to = c(paste0("C",seq(0,3))))
scobj.query.batch$celltype <- factor(scobj.query.batch$celltype,levels=c(paste0("C",seq(0,3))))
saveRDS(scobj.query.batch, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.query.scobj.SOT.V4.Rds")
DimPlot(scobj.query.batch, reduction="umap", label = TRUE,pt.size=0.5, group.by="celltype")

## merge two datasets
scobj.batch <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")
scobj.query.batch <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.query.scobj.SOT.V4.Rds")
scobj.combined <- merge(x=scobj.batch, y=scobj.query.batch, add.cell.ids = c("Seurat1", "Seurat2"))
scobj.combined <- NormalizeData(scobj.combined, normalization.method = "LogNormalize", scale.factor = 10000)
scobj.combined <- FindVariableFeatures(scobj.combined, selection.method = "vst", verbose = T, nfeatures = 3000)
scobj.combined <- CellCycleScoring(scobj.combined, s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)
scobj.combined <- ScaleData(scobj.combined, features = rownames(scobj.combined))

scobj.combined <- RunPCA(scobj.combined, features=hvg.genes.add, verbose = F, npcs = 100)
ElbowPlot(scobj.combined)
ElbowPlot(scobj.combined, reduction = "pca",ndims=15)
DimPlot(scobj.combined, reduction = "pca",group.by="orig.ident")
scobj.combined <- RunUMAP(scobj.combined, reduction = "pca",dims=1:15)
DimPlot(scobj.combined, reduction = "umap",group.by="orig.ident")
#batch
scobj.combined.batch <- RunHarmony(scobj.combined, group.by.vars="orig.ident")
scobj.combined.batch <- RunUMAP(scobj.combined.batch, reduction = "harmony",dims=1:15)
DimPlot(scobj.combined.batch, reduction="umap",group.by="orig.ident",pt.size=0.5)
#cluster
DimPlot(scobj.combined.batch, reduction="umap", label = TRUE,pt.size=0.5, group.by="celltype")


#====== 2 ReName celltype and plot ============================
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")
#sce$celltype <- plyr::mapvalues(x=sce$population,
#                                from = c("0","6","5","3","7","4","1","2","8"),
#                                to = c(paste0("C",seq(0,8))))
#sce$celltype <- factor(sce$celltype,levels=c(paste0("C",seq(0,8))))
#saveRDS(sce, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")

cell <- sce@meta.data %>% mutate(BC=rownames(.)) %>% select(BC, orig.ident, celltype) 
cell$BC <- gsub("-1","",cell$BC)
cell$BC <- gsub(".*SC_", "BC_", cell$BC)
cell$orig.ident <- gsub("10X_","",cell$orig.ident)
cell$orig.ident <- gsub("H9_","hESC_",cell$orig.ident)
unique(cell$orig.ident)
table(cell$orig.ident)
#write.csv(cell, "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.all.cell.infos.csv", row.names=F)

#==plot===========================================================================================
#Figure 4B UMAP sample
sample <- c("10X_H9_A7","10X_H9_C9","10X_H9_C12","10X_H9_D3","10X_HSC_B7","10X_HSC_E12","10X_HSC_F12")
color <- c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#FB9A99","#6A3D9A","#FDBF6F")
color <- setNames(color, sample)

sce$orig.ident <- factor(sce$orig.ident, levels = sample)
DimPlot(sce,reduction="umap",group.by="orig.ident",pt.size=0.1,cols = color) +
  labs(x="UMAP1",y="UMAP2",title="")


#Figure 4C UMAP cluster
clusterCols <- c("#B2DF8A","#3DA42D","#FE9999","#A4CEDE","#FDBE6F","#FF9116","#BBBFE0","#673A95","#B25A2C")
DimPlot(sce, reduction="umap",label = T, cols=clusterCols, pt.size=0.1,group.by="celltype")+labs(x="UMAP1",y="UMAP2",title="")

#Figure 4D Cluster distribution per sample
cluster <- sce@meta.data[,c(1,19)]
cluster$orig.ident <- plyr::mapvalues(cluster$orig.ident,
                                      from=c("10X_H9_A7","10X_H9_C9","10X_H9_C12","10X_H9_D3","10X_HSC_B7","10X_HSC_E12","10X_HSC_F12"),
                                      to=c("A7_hESC","C9_hESC","C12_hESC","D3_hESC","B7_HSC","E12_HSC","F12_HSC"))
res1 <- cluster %>% group_by(orig.ident) %>% dplyr::count(celltype) %>% ungroup()
res2 <- res1 %>% dplyr::group_by(orig.ident) %>% dplyr::summarise(n_sum = sum(n)) %>% ungroup()
res <- merge(res1,res2,by='orig.ident') %>% mutate(ratio=n/n_sum) %>% select(1,2,5)
res$orig.ident <- factor(res$orig.ident, levels=c("A7_hESC","C9_hESC","C12_hESC","D3_hESC","B7_HSC","E12_HSC","F12_HSC"))

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4/SOT/numBalance.hsc.hesc.cluster.distribution.pdf",width=4,height=5)
clusterCols <- c("#B2DF8A","#3DA42D","#FE9999","#A4CEDE","#FDBE6F","#FF9116","#BBBFE0","#673A95","#B25A2C")
ggplot(res,aes(x=orig.ident,y=ratio,fill=celltype))+geom_bar(stat="identity", width=0.6, position="fill",alpha=0.8)+
  labs(x="",y="")+scale_fill_manual(values=clusterCols)+ #theme_classic()+ 
  #scale_y_continuous(expand=c(0,0))+
  theme(panel.border=element_rect(fill='transparent',color='black'),panel.background=element_blank(),
        axis.ticks=element_line(color="black"),axis.ticks.length=unit(0.13,"cm"),
        axis.title.y = element_text(size=12,color="black"),axis.text.x=element_text(size=10,color='black',angle=30,hjust=1,vjust=1),
        axis.text.y=element_text(size=10,color="black"))

#Figure 4E DEGs exp ration per cluster  (C1-C10)
p.sc <- DotPlot(sce,group.by="celltype",features=c("D4_sc_Score","D6_sc_Score","D8_sc_Score","D12_sc_Score"),
                cols = "RdYlBu",dot.scale=7)+labs(x="",y="")+ theme(axis.text.x=element_text(angle=30,hjust=1),legend.position="none")
p.array <- DotPlot(sce,group.by="celltype",features=c("D0_array_Score","D12_array_Score"),
                   cols = "RdYlBu",dot.scale=7)+labs(x="",y="")+ theme(axis.text.x=element_text(angle=30,hjust=1))
p.sc+p.array


#Figure 4F marker gene expression per cluster
marker_genelist.new <- c("POU5F1","NANOG","SOX2","ESRG","GDF3","THY1",
                         "GATA6","RGS5","GPRC5C","FABP5","CAP2","HOXB2",
                         "ICAM1","FGF1",  "CCL5","IGFBP6","TIMP1","SDC4",
                         "PDLIM4","FSTL3","TNFRSF12A","TIMP2","EMP3","PDGFA",
                         "CTHRC1","COL4A1","PDGFB","IGFBP3","MYL9","FBLN1","CRABP2","TAGLN","FN1",
                         "COL1A1","LUM","ACTA2","CCL2","COL3A1","TGFBR1","GPC3",
                         "GDF15","HOPX","CD55","TNFRSF1B","HMOX1","CDKN1C")
DotPlot(sce, features=marker_genelist.new, group.by="celltype",cols ="RdYlBu",dot.scale=7)+ labs(x="",y="")+
  theme(axis.text.x=element_text(angle=30,hjust=1,size=10),axis.text.y=element_text(size=10),
        legend.title = element_text(size=10),legend.text = element_text(size=10),legend.position = "right")


