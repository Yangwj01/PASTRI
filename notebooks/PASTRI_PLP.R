rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  #library(tidyverse)
  library(PASTRI)
  library(parallel)
  library(ape)
  library(ggtree)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
})

## ==================== source and define functions ===========================================
plot.TreeAndCelltype <- function(tree, all.info.df, sample.name, color_value, alphas){
  # 1. extract all node celltype infos of all nodes
  node.celltype.info <- all.info.df %>% filter(sample == sample.name & celltype != "inter")
  # 2. note all multi cell nodes
  select.node.celltype <- node.celltype.info %>% ungroup() %>%
    mutate(celltype=ifelse(cellNum>1, "multi_nodes", celltype)) %>%
    select(label, sample, cellNum, celltype, type) %>% unique() %>% as.data.frame()
  rownames(select.node.celltype) <- select.node.celltype$label
  
  # 4. set use colors
  celltype.colors <- structure(c("#1F78B4","#B2DF8A","#A6CEE3","#33A02C","#E31A1C","#FF7F00","#FDBF6F","#B15928","#6A3D9A","#CAB2D6","#FFFF99","#FB9A99","#666666"),
                               names=as.character(c(paste0("C",seq(1,10)), "R1", "R2", "multi_nodes")))
  # 5. only keep the celltype cols
  select.node.celltype.use <- as.data.frame(select.node.celltype[, "celltype", drop=F])
  rownames(select.node.celltype.use) <- select.node.celltype$label
  select.node.celltype.use$celltype <- as.factor(select.node.celltype.use$celltype)
  
  # plot tree and celltypes
  p <- ggtree(tree,branch.length="none",layout="circular",size=0.1,color=color_value,alpha=alphas)+xlim(-10, NA) #circular or radial
  p <- rotate_tree(p, 30)
  p <- 
    gheatmap(p, select.node.celltype.use, offset=0.05, width=0.1, colnames_angle=120, colnames_offset_y = .25, colnames = FALSE,color=NA) + xlim(-2, NA) +
    scale_fill_manual(name="cell types",values = celltype.colors)
  return(p)
}


## ============ Run ===========================================================================
#== 0 get cell infos
cell_info <- readRDS("/mnt/data4/disk/zhxyu8/lung_progenitor_diff_lineage/Data_upload/all_cbrad5_GS_hesc_tree_dataframe_modify_new.Rds")
cell_info <- cell_info %>% filter(type=="leaf")
cell_info$celltype <- plyr::mapvalues(x=cell_info$celltype,
                                      from = c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                      to = c(paste0("C",seq(1,10)), "R1", "R2"))
cell_infos <- mclapply(seq(length(unique(cell_info$sample))),function(s){
  sample <- unique(cell_info$sample)
  sel.sample <- sample[s]
  cell.info <- cell_info %>% filter(sample==sel.sample) %>% group_by(label) %>% add_count(label,name="cellNum") %>% ungroup()
  return(cell.info)
}) %>% bind_rows()

all.cell_infos <- cell_infos %>% select(label,BC,sample,celltype,cellNum)
write.csv(all.cell_infos,"/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/all.cell_infos.csv",row.names=F)

##Figure 6A plot example tree
A1.tree <- read.tree("/mnt/data4/disk/zhxyu8/lung_progenitor_diff_lineage/Data_upload/A1-CBRAD5.nwk")
p.all.A1 <- plot.TreeAndCelltype(tree=A1.tree, all.info.df=cell_infos, sample.name="A1-CBRAD5", 
                                 color_value="black", alphas=0.3)
pdf("/mnt/data5/disk/yangwj/Result_plots/Figure6/Figure6A.pdf",width=6,height=6) 
p.all.A1
dev.off()


#== 1 get optimal transition
#calculate_lca_depths
tree_path <- "/mnt/data4/disk/zhxyu8/lung_progenitor_diff_lineage/Data_upload/"
A1.node.depth <- calculate_lca_depths(paste0(tree_path, "A1-CBRAD5.nwk"))
G2.node.depth <- calculate_lca_depths(paste0(tree_path, "G2-CBRAD5.nwk"))
G11.node.depth <- calculate_lca_depths(paste0(tree_path, "G11-CBRAD5.nwk"))

save(A1.node.depth, G2.node.depth, G11.node.depth,
     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/PLP.Node.depth.Rda")

#load some files
cell <- read.csv("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/all.cell_infos.csv")
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/PLP.Node.depth.Rda")

#get optimal trans
Sample <- c("A1", "G2", "G11") 
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

df_mrca.mheight <- lapply(seq(length(Sample)), function(S){
  #S <- 3
  sel.Sample <- Sample[S]
  print(sel.Sample)
  sel.Onenode.depth <- paste0(sel.Sample,".node.depth") %>% get()
  sel.sample.name <- paste0(sel.Sample, "-CBRAD5")
  
  ##get cell info
  sel.cell <- cell %>% filter(sample==sel.sample.name) %>% select(label, BC, celltype, cellNum)
  colnames(sel.cell) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
  #sel.cell %>% filter(cellNum==1) %>% pull(celltype) %>% table()
  ##
  cell_type <- c("C5","C6","C7","C8","C9","C10")
  Bound_matrix <- matrix(c(Inf,Inf,Inf,Inf,Inf,Inf, 
                           0,Inf,Inf,Inf,Inf,Inf, 
                           0,0,Inf,Inf,Inf,Inf, 
                           0,0,0,Inf,Inf,Inf, 
                           0,0,0,0,Inf,0, 
                           0,0,0,0,0,Inf), nrow=length(cell_type), byrow=FALSE)
  
  ##Infer optimal constrained transition matrix 
  df_depth <- lapply(seq(length(depth)), function(d){
    #d <- 1
    sel.depth <- depth[d]
    print(sel.depth)
    
    # Call get_optimal_transition_matrix() with the specified parameters
    mrca.mheight_optimal_results_d <- get_optimal_transition_matrix(
      node_pair_depth = sel.Onenode.depth,
      cell_info = sel.cell,
      Sel_u = "lca_normalized_height",
      fi_depth = sel.depth,
      Bound_Matrix = Bound_matrix,
      cell_type_list = cell_type,
      mc.cores = 80
    )
    mrca.mheight_optimal_results_d1 <- mrca.mheight_optimal_results_d$optimal_norm_df_dataframe %>%
      mutate(sample=sel.sample.name)
    return(mrca.mheight_optimal_results_d1)
  }) %>% bind_rows()
  ##
  out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/"
  saveRDS(df_depth, file= paste0(out_path, sel.Sample,".op.trans",".Rds"))
  print("Successful")
  return(df_depth)
})

#== 2 compare correlation between samples
#rm(list=ls())
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
A1.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/A1.op.trans.Rds")
G2.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/G2.op.trans.Rds")
G11.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/G11.op.trans.Rds")

cor.res <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  sel.depth.1 <- paste0("d_", sel.depth)
  print(sel.depth.1)
  ##
  irreversible.type <- c("C5_C5","C5_C6","C5_C7","C5_C8","C5_C9","C5_C10",
                         "C6_C6","C6_C7","C6_C8","C6_C9","C6_C10",
                         "C7_C7","C7_C8","C7_C9","C7_C10",
                         "C8_C8","C8_C9","C8_C10",
                         "C9_C9", "C10_C10")
  sel.A1.op <- A1.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>%
    select(type_combined_new, start_cell, end_cell, norm_optimal_T)
  colnames(sel.A1.op) <- c("type.combined.new","start.cell","end.cell","A1.norm.optimal.T")
  sel.G2.op <- G2.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>% 
    select(type_combined_new, norm_optimal_T)
  colnames(sel.G2.op) <- c("type.combined.new","G2.norm.optimal.T")
  sel.G11.op <- G11.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>%
    select(type_combined_new, norm_optimal_T)
  colnames(sel.G11.op) <- c("type.combined.new","G11.norm.optimal.T")
  all <- left_join(left_join(sel.A1.op, sel.G2.op, by="type.combined.new"), sel.G11.op, by="type.combined.new")
  ##compare
  A1_G2.cor <- cor.test(all$A1.norm.optimal.T, all$G2.norm.optimal.T, method = "pearson")
  A1_G2.rho <- cor.test(all$A1.norm.optimal.T, all$G2.norm.optimal.T, method = "spearman")
  A1_G2.fdist <- dist(rbind(all$A1.norm.optimal.T, all$G2.norm.optimal.T), method = "euclidean") %>% as.vector()
  A1_G2.mhdist <- dist(rbind(all$A1.norm.optimal.T, all$G2.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  A1_G11.cor <- cor.test(all$A1.norm.optimal.T, all$G11.norm.optimal.T, method = "pearson")
  A1_G11.rho <- cor.test(all$A1.norm.optimal.T, all$G11.norm.optimal.T, method = "spearman")
  A1_G11.fdist <- dist(rbind(all$A1.norm.optimal.T, all$G11.norm.optimal.T), method = "euclidean") %>% as.vector()
  A1_G11.mhdist <- dist(rbind(all$A1.norm.optimal.T, all$G11.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  G2_G11.cor <- cor.test(all$G2.norm.optimal.T, all$G11.norm.optimal.T, method = "pearson")
  G2_G11.rho <- cor.test(all$G2.norm.optimal.T, all$G11.norm.optimal.T, method = "spearman")
  G2_G11.fdist <- dist(rbind(all$G2.norm.optimal.T, all$G11.norm.optimal.T), method = "euclidean") %>% as.vector()
  G2_G11.mhdist <- dist(rbind(all$G2.norm.optimal.T, all$G11.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  ##
  return(data.frame(depth=sel.depth,
                    times="real", 
                    stringsAsFactors =F,
                    A1_G2.cor=A1_G2.cor$estimate, A1_G2.cor.pval=A1_G2.cor$p.value,
                    A1_G2.rho=A1_G2.rho$estimate, A1_G2.rho.pval=A1_G2.rho$p.value,
                    A1_G2.fdist=A1_G2.fdist, A1_G2.mhdist=A1_G2.mhdist,
                    
                    A1_G11.cor=A1_G11.cor$estimate, A1_G11.cor.pval=A1_G11.cor$p.value,
                    A1_G11.rho=A1_G11.rho$estimate, A1_G11.rho.pval=A1_G11.rho$p.value,
                    A1_G11.fdist=A1_G11.fdist, A1_G11.mhdist=A1_G11.mhdist,
                    
                    G2_G11.cor=G2_G11.cor$estimate, G2_G11.cor.pval=G2_G11.cor$p.value,
                    G2_G11.rho=G2_G11.rho$estimate, G2_G11.rho.pval=G2_G11.rho$p.value,
                    G2_G11.fdist=G2_G11.fdist, G2_G11.mhdist=G2_G11.mhdist))
}) %>% bind_rows() %>%
  mutate(mean.cor=(A1_G2.cor+A1_G11.cor+G2_G11.cor)/3,
         mean.rho=(A1_G2.rho+A1_G11.rho+G2_G11.rho)/3,
         mean.fdist=(A1_G2.fdist+A1_G11.fdist+G2_G11.fdist)/3,
         mean.mhdist=(A1_G2.mhdist+A1_G11.mhdist+G2_G11.mhdist)/3)
cor.res.new <- lapply(seq(length(depth)), function(d){
  #d <- 6
  sel.depth <- depth[d]
  print(sel.depth)
  
  sel.cor.res <- cor.res %>% filter(depth==sel.depth) %>%
    mutate(mean.cor=(A1_G2.cor+A1_G11.cor+G2_G11.cor)/3,
           mean.fdist=(A1_G2.fdist+A1_G11.fdist+G2_G11.fdist)/3,
           sd.fdist=sd(c(A1_G2.fdist,A1_G11.fdist,G2_G11.fdist)), se.fdist=sd.fdist/sqrt(3))
  return(sel.cor.res)
}) %>% bind_rows() %>% mutate(diff1=mean.fdist-se.fdist, diff2=mean.fdist+se.fdist)

saveRDS(cor.res, "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/compare_real.Rds")

## Figure 6C
p.cor.res <- lapply(seq(length(depth)), function(d){
  #d <- 5
  sel.depth <- depth[d]
  sel.depth.1 <- paste0("d_", sel.depth)
  print(sel.depth.1)
  ##
  irreversible.type <- c("C5_C5","C5_C6","C5_C7","C5_C8","C5_C9","C5_C10",
                         "C6_C6","C6_C7","C6_C8","C6_C9","C6_C10",
                         "C7_C7","C7_C8","C7_C9","C7_C10",
                         "C8_C8","C8_C9","C8_C10",
                         "C9_C9","C10_C10")
  sel.A1.op <- A1.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>%
    select(type_combined_new, start_cell, end_cell, norm_optimal_T)
  colnames(sel.A1.op) <- c("type.combined.new","start.cell","end.cell","A1.norm.optimal.T")
  sel.G2.op <- G2.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>% 
    select(type_combined_new, norm_optimal_T)
  colnames(sel.G2.op) <- c("type.combined.new","G2.norm.optimal.T")
  sel.G11.op <- G11.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>% 
    select(type_combined_new, norm_optimal_T)
  colnames(sel.G11.op) <- c("type.combined.new","G11.norm.optimal.T")
  all <- left_join(left_join(sel.A1.op, sel.G2.op, by="type.combined.new"), sel.G11.op, by="type.combined.new")
  all$start.cell <- factor(all$start.cell, levels = c("C5","C6","C7","C8","C9","C10"))
  all$end.cell <- factor(all$end.cell, levels = c("C5","C6","C7","C8","C9","C10"))
  ##
  sel.cor.res <- cor.res %>% filter(depth==sel.depth)
  
  ##
  p.A1_G2 <- ggplot(all, aes(x = A1.norm.optimal.T, y = G2.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=3)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$A1_G2.cor, 2), " P = ", log10(sel.cor.res$A1_G2.cor.pval), 
                      "\neuclidean distance = ", round(sel.cor.res$A1_G2.fdist, 2)),
         x="PASTRI inferred transition rate of A1",y="PASTRI inferred transition rate of G2")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+
    xlim(0,1)+ylim(0,1)+scale_shape_manual(values=c(0,6,7,15,16,17))+
    scale_color_manual(values=c("#E31A1C","#FF7F00","#FDBF6F","#B15928","#6A3D9A","#CAB2D6"))+
    theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
          legend.position = "none")
  p.A1_G11 <- ggplot(all, aes(x = A1.norm.optimal.T, y = G11.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=3)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$A1_G11.cor, 2), " P = ", log10(sel.cor.res$A1_G11.cor.pval), 
                      "\neuclidean distance = ", round(sel.cor.res$A1_G11.fdist, 2)),
         x="PASTRI inferred transition rate of A1",y="PASTRI inferred transition rate of G11")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+
    xlim(0,1)+ylim(0,1)+scale_shape_manual(values=c(0,6,7,15,16,17))+
    scale_color_manual(values=c("#E31A1C","#FF7F00","#FDBF6F","#B15928","#6A3D9A","#CAB2D6"))+
    theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
          legend.position = "none")
  p.G2_G11 <- ggplot(all, aes(x = G2.norm.optimal.T, y = G11.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=3)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$G2_G11.cor, 2), " P = ", log10(sel.cor.res$G2_G11.cor.pval),
                      "\neuclidean distance = ", round(sel.cor.res$G2_G11.fdist, 2)),
         x="PASTRI inferred transition rate of G2",y="PASTRI inferred transition rate of G11")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+
    xlim(0,1)+ylim(0,1)+scale_shape_manual(values=c(0,6,7,15,16,17))+
    scale_color_manual(values=c("#E31A1C","#FF7F00","#FDBF6F","#B15928","#6A3D9A","#CAB2D6"))+
    theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
          legend.text=element_text(size=10,color="black"))
  p <- p.A1_G2+p.A1_G11+p.G2_G11
  return(p)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure6/PASTRI/Figure6C.pdf",width=12,height=5) 
p.cor.res[[5]]
dev.off()


#== 3 PASTRI accuracy
## 3.1 shuffle cellInfo and re-calculate transition rate by PASTRI
rm(list=ls())
cell <- read.csv("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/all.cell_infos.csv")
colnames(cell) <- c("nodeLabel", "cell_ID", "sample", "celltype", "cellNum")
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/PLP.Node.depth.Rda")
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
sel.sample <- "A1"
other.sample.1 <- "G2"
other.sample.2 <- "G11"
n <- 1000

df_shuffle <- lapply(seq(n), function(i){
  #i <- 1
  repeat{
    cat("times", i, "\n")
    sel.Onenode.depth <- paste0(sel.sample,".node.depth") %>% get()
    other1.Onenode.depth <- paste0(other.sample.1,".node.depth") %>% get()
    other2.Onenode.depth <- paste0(other.sample.2,".node.depth") %>% get()
    ##load sel sample celltype infos and random tree leaf
    sel.celltype <- subset(cell, cell$sample==paste0(sel.sample, "-CBRAD5"))
    other1.celltype <- subset(cell, cell$sample==paste0(other.sample.1, "-CBRAD5"))
    other2.celltype <- subset(cell, cell$sample==paste0(other.sample.2, "-CBRAD5"))
    ##random tree leaf
    sel.sim.One.celltype <- sel.celltype %>% filter(cellNum==1) %>% 
      #mutate(celltype=sample(celltype, size=nrow(.), replace=F))
      mutate(celltype=sample(sel.celltype$celltype, size=nrow(.), replace=F))
    #table(sel.sim.One.celltype$celltype)
    other1.sim.One.celltype <- other1.celltype %>% filter(cellNum==1) %>% 
      #mutate(celltype=sample(celltype, size=nrow(.), replace=F))
      mutate(celltype=sample(other1.celltype$celltype, size=nrow(.), replace=F))
    other2.sim.One.celltype <- other2.celltype %>% filter(cellNum==1) %>% 
      #mutate(celltype=sample(celltype, size=nrow(.), replace=F))
      mutate(celltype=sample(other2.celltype$celltype, size=nrow(.), replace=F))
    
    #
    cell_type <- c("C5","C6","C7","C8","C9","C10")
    Bound_matrix <- matrix(c(Inf,Inf,Inf,Inf,Inf,Inf, 
                             0,Inf,Inf,Inf,Inf,Inf, 
                             0,0,Inf,Inf,Inf,Inf, 
                             0,0,0,Inf,Inf,Inf, 
                             0,0,0,0,Inf,0, 
                             0,0,0,0,0,Inf), nrow=length(cell_type), byrow=FALSE)
    
    # Check if conditions are met
    df_depth_try <- tryCatch({  
      lapply(seq(length(depth)), function(d){
        #d <- 2
        sel.depth <- depth[d]
        print(sel.depth)
        ##random
        sel_optimal_results_d <- get_optimal_transition_matrix(node_pair_depth = sel.Onenode.depth,
                                                               cell_info = sel.sim.One.celltype,
                                                               Sel_u = "lca_normalized_height",
                                                               fi_depth = sel.depth,
                                                               Bound_Matrix = Bound_matrix,
                                                               cell_type_list = cell_type,
                                                               mc.cores = 80) %>% .$optimal_norm_df_dataframe %>% mutate(sample=sel.sample)
        other_1_optimal_results_d <- get_optimal_transition_matrix(node_pair_depth = other1.Onenode.depth,
                                                                   cell_info = other1.sim.One.celltype,
                                                                   Sel_u = "lca_normalized_height",
                                                                   fi_depth = sel.depth,
                                                                   Bound_Matrix = Bound_matrix,
                                                                   cell_type_list = cell_type,
                                                                   mc.cores = 80) %>% .$optimal_norm_df_dataframe %>% mutate(sample=other.sample.1)
        other_2_optimal_results_d <- get_optimal_transition_matrix(node_pair_depth = other2.Onenode.depth,
                                                                   cell_info = other2.sim.One.celltype,
                                                                   Sel_u = "lca_normalized_height",
                                                                   fi_depth = sel.depth,
                                                                   Bound_Matrix = Bound_matrix,
                                                                   cell_type_list = cell_type,
                                                                   mc.cores = 80) %>% .$optimal_norm_df_dataframe %>% mutate(sample=other.sample.2)
        ##real
        sel_real_A1 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/A1.op.trans.Rds") %>%
          filter(depth==paste0("d_",sel.depth))
        sel_real_G2 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/G2.op.trans.Rds") %>%
          filter(depth==paste0("d_",sel.depth))
        sel_real_G11 <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/G11.op.trans.Rds") %>%
          filter(depth==paste0("d_",sel.depth))
        
        ##compare correlation between samples
        ##
        irreversible.type <- c("C5_C5","C5_C6","C5_C7","C5_C8","C5_C9","C5_C10",
                               "C6_C6","C6_C7","C6_C8","C6_C9","C6_C10",
                               "C7_C7","C7_C8","C7_C9","C7_C10",
                               "C8_C8","C8_C9","C8_C10",
                               "C9_C9", "C10_C10")
        #random
        sel.op <- sel_optimal_results_d %>% filter(type_combined_new %in% irreversible.type) %>%
          select(type_combined_new, start_cell, end_cell, norm_optimal_T)
        colnames(sel.op) <- c("type.combined.new","start.cell","end.cell", paste0(sel.sample, ".norm.optimal.T"))
        other_1.op <- other_1_optimal_results_d %>% filter(type_combined_new %in% irreversible.type) %>% 
          select(type_combined_new, norm_optimal_T)
        colnames(other_1.op) <- c("type.combined.new", paste0(other.sample.1, ".norm.optimal.T"))
        other_2.op <- other_2_optimal_results_d %>% filter(type_combined_new %in% irreversible.type) %>% 
          select(type_combined_new, norm_optimal_T)
        colnames(other_2.op) <- c("type.combined.new", paste0(other.sample.2, ".norm.optimal.T"))
        #real
        A1.op <- sel_real_A1 %>% filter(type_combined_new %in% irreversible.type) %>% select(type_combined_new, norm_optimal_T)
        colnames(A1.op) <- c("type.combined.new", "real.A1.norm.optimal.T")
        G2.op <- sel_real_G2 %>% filter(type_combined_new %in% irreversible.type) %>% select(type_combined_new, norm_optimal_T)
        colnames(G2.op) <- c("type.combined.new", "real.G2.norm.optimal.T")
        G11.op <- sel_real_G11 %>% filter(type_combined_new %in% irreversible.type) %>% select(type_combined_new, norm_optimal_T)
        colnames(G11.op) <- c("type.combined.new", "real.G11.norm.optimal.T")
        
        all <- full_join(full_join(full_join(sel.op, other_1.op, by="type.combined.new"),
                                   full_join(other_2.op, A1.op, by="type.combined.new"), by="type.combined.new"),
                         full_join(G2.op, G11.op, by="type.combined.new"), by="type.combined.new")
        ##compare
        sim_A1_real_G2.cor <- cor.test(all$A1.norm.optimal.T, all$real.G2.norm.optimal.T, method = "pearson")
        sim_A1_real_G2.rho <- cor.test(all$A1.norm.optimal.T, all$real.G2.norm.optimal.T, method = "spearman")
        sim_A1_real_G2.fdist <- dist(rbind(all$A1.norm.optimal.T, all$real.G2.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_A1_real_G2.mhdist <- dist(rbind(all$A1.norm.optimal.T, all$real.G2.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_A1_real_G11.cor <- cor.test(all$A1.norm.optimal.T, all$real.G11.norm.optimal.T, method = "pearson")
        sim_A1_real_G11.rho <- cor.test(all$A1.norm.optimal.T, all$real.G11.norm.optimal.T, method = "spearman")
        sim_A1_real_G11.fdist <- dist(rbind(all$A1.norm.optimal.T, all$real.G11.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_A1_real_G11.mhdist <- dist(rbind(all$A1.norm.optimal.T, all$real.G11.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_G2_real_A1.cor <- cor.test(all$G2.norm.optimal.T, all$real.A1.norm.optimal.T, method = "pearson")
        sim_G2_real_A1.rho <- cor.test(all$G2.norm.optimal.T, all$real.A1.norm.optimal.T, method = "spearman")
        sim_G2_real_A1.fdist <- dist(rbind(all$G2.norm.optimal.T, all$real.A1.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_G2_real_A1.mhdist <- dist(rbind(all$G2.norm.optimal.T, all$real.A1.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_G2_real_G11.cor <- cor.test(all$G2.norm.optimal.T, all$real.G11.norm.optimal.T, method = "pearson")
        sim_G2_real_G11.rho <- cor.test(all$G2.norm.optimal.T, all$real.G11.norm.optimal.T, method = "spearman")
        sim_G2_real_G11.fdist <- dist(rbind(all$G2.norm.optimal.T, all$real.G11.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_G2_real_G11.mhdist <- dist(rbind(all$G2.norm.optimal.T, all$real.G11.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_G11_real_A1.cor <- cor.test(all$G11.norm.optimal.T, all$real.A1.norm.optimal.T, method = "pearson")
        sim_G11_real_A1.rho <- cor.test(all$G11.norm.optimal.T, all$real.A1.norm.optimal.T, method = "spearman")
        sim_G11_real_A1.fdist <- dist(rbind(all$G11.norm.optimal.T, all$real.A1.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_G11_real_A1.mhdist <- dist(rbind(all$G11.norm.optimal.T, all$real.A1.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_G11_real_G2.cor <- cor.test(all$G11.norm.optimal.T, all$real.G2.norm.optimal.T, method = "pearson")
        sim_G11_real_G2.rho <- cor.test(all$G11.norm.optimal.T, all$real.G2.norm.optimal.T, method = "spearman")
        sim_G11_real_G2.fdist <- dist(rbind(all$G11.norm.optimal.T, all$real.G2.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_G11_real_G2.mhdist <- dist(rbind(all$G11.norm.optimal.T, all$real.G2.norm.optimal.T), method = "manhattan") %>% as.vector()
        ##
        sim_A1_G2.cor <- cor.test(all$A1.norm.optimal.T, all$G2.norm.optimal.T, method = "pearson")
        sim_A1_G2.rho <- cor.test(all$A1.norm.optimal.T, all$G2.norm.optimal.T, method = "spearman")
        sim_A1_G2.fdist <- dist(rbind(all$A1.norm.optimal.T, all$G2.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_A1_G2.mhdist <- dist(rbind(all$A1.norm.optimal.T, all$G2.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_A1_G11.cor <- cor.test(all$A1.norm.optimal.T, all$G11.norm.optimal.T, method = "pearson")
        sim_A1_G11.rho <- cor.test(all$A1.norm.optimal.T, all$G11.norm.optimal.T, method = "spearman")
        sim_A1_G11.fdist <- dist(rbind(all$A1.norm.optimal.T, all$G11.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_A1_G11.mhdist <- dist(rbind(all$A1.norm.optimal.T, all$G11.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_G2_G11.cor <- cor.test(all$G2.norm.optimal.T, all$G11.norm.optimal.T, method = "pearson")
        sim_G2_G11.rho <- cor.test(all$G2.norm.optimal.T, all$G11.norm.optimal.T, method = "spearman")
        sim_G2_G11.fdist <- dist(rbind(all$G2.norm.optimal.T, all$G11.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_G2_G11.mhdist <- dist(rbind(all$G2.norm.optimal.T, all$G11.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        ##
        return(data.frame(depth=sel.depth,
                          times=i, 
                          stringsAsFactors =F,
                          sim_A1_real_G2.cor=sim_A1_real_G2.cor$estimate, sim_A1_real_G2.cor.pval=sim_A1_real_G2.cor$p.value,
                          sim_A1_real_G2.rho=sim_A1_real_G2.rho$estimate, sim_A1_real_G2.rho.pval=sim_A1_real_G2.rho$p.value,
                          sim_A1_real_G2.fdist=sim_A1_real_G2.fdist, sim_A1_real_G2.mhdist=sim_A1_real_G2.mhdist,
                          
                          sim_A1_real_G11.cor=sim_A1_real_G11.cor$estimate, sim_A1_real_G11.cor.pval=sim_A1_real_G11.cor$p.value,
                          sim_A1_real_G11.rho=sim_A1_real_G11.rho$estimate, sim_A1_real_G11.rho.pval=sim_A1_real_G11.rho$p.value,
                          sim_A1_real_G11.fdist=sim_A1_real_G11.fdist, sim_A1_real_G11.mhdist=sim_A1_real_G11.mhdist,
                          
                          sim_G2_real_A1.cor=sim_G2_real_A1.cor$estimate, sim_G2_real_A1.cor.pval=sim_G2_real_A1.cor$p.value,
                          sim_G2_real_A1.rho=sim_G2_real_A1.rho$estimate, sim_G2_real_A1.rho.pval=sim_G2_real_A1.rho$p.value,
                          sim_G2_real_A1.fdist=sim_G2_real_A1.fdist, sim_G2_real_A1.mhdist=sim_G2_real_A1.mhdist,
                          
                          sim_G2_real_G11.cor=sim_G2_real_G11.cor$estimate, sim_G2_real_G11.cor.pval=sim_G2_real_G11.cor$p.value,
                          sim_G2_real_G11.rho=sim_G2_real_G11.rho$estimate, sim_G2_real_G11.rho.pval=sim_G2_real_G11.rho$p.value,
                          sim_G2_real_G11.fdist=sim_G2_real_G11.fdist, sim_G2_real_G11.mhdist=sim_G2_real_G11.mhdist,
                          
                          sim_G11_real_A1.cor=sim_G11_real_A1.cor$estimate, sim_G11_real_A1.cor.pval=sim_G11_real_A1.cor$p.value,
                          sim_G11_real_A1.rho=sim_G11_real_A1.rho$estimate, sim_G11_real_A1.rho.pval=sim_G11_real_A1.rho$p.value,
                          sim_G11_real_A1.fdist=sim_G11_real_A1.fdist, sim_G11_real_A1.mhdist=sim_G11_real_A1.mhdist,
                          
                          sim_G11_real_G2.cor=sim_G11_real_G2.cor$estimate, sim_G11_real_G2.cor.pval=sim_G11_real_G2.cor$p.value,
                          sim_G11_real_G2.rho=sim_G11_real_G2.rho$estimate, sim_G11_real_G2.rho.pval=sim_G11_real_G2.rho$p.value,
                          sim_G11_real_G2.fdist=sim_G11_real_G2.fdist, sim_G11_real_G2.mhdist=sim_G11_real_G2.mhdist,
                          #
                          sim_A1_G2.cor=sim_A1_G2.cor$estimate, sim_A1_G2.cor.pval=sim_A1_G2.cor$p.value,
                          sim_A1_G2.rho=sim_A1_G2.rho$estimate, sim_A1_G2.rho.pval=sim_A1_G2.rho$p.value,
                          sim_A1_G2.fdist=sim_A1_G2.fdist, sim_A1_G2.mhdist=sim_A1_G2.mhdist,
                          
                          sim_A1_G11.cor=sim_A1_G11.cor$estimate, sim_A1_G11.cor.pval=sim_A1_G11.cor$p.value,
                          sim_A1_G11.rho=sim_A1_G11.rho$estimate, sim_A1_G11.rho.pval=sim_A1_G11.rho$p.value,
                          sim_A1_G11.fdist=sim_A1_G11.fdist, sim_A1_G11.mhdist=sim_A1_G11.mhdist,
                          
                          sim_G2_G11.cor=sim_G2_G11.cor$estimate, sim_G2_G11.cor.pval=sim_G2_G11.cor$p.value,
                          sim_G2_G11.rho=sim_G2_G11.rho$estimate, sim_G2_G11.rho.pval=sim_G2_G11.rho$p.value,
                          sim_G2_G11.fdist=sim_G2_G11.fdist, sim_G2_G11.mhdist=sim_G2_G11.mhdist,
                          ##
                          mean.cor=mean(c(sim_A1_real_G2.cor$estimate,sim_A1_real_G11.cor$estimate,
                                          sim_G2_real_A1.cor$estimate,sim_G2_real_G11.cor$estimate,
                                          sim_G11_real_A1.cor$estimate,sim_G11_real_G2.cor$estimate)),
                          mean.rho=mean(c(sim_A1_real_G2.rho$estimate,sim_A1_real_G11.rho$estimate,
                                          sim_G2_real_A1.rho$estimate,sim_G2_real_G11.rho$estimate,
                                          sim_G11_real_A1.rho$estimate,sim_G11_real_G2.rho$estimate)),
                          mean.fdist=mean(c(sim_A1_real_G2.fdist,sim_A1_real_G11.fdist,
                                            sim_G2_real_A1.fdist,sim_G2_real_G11.fdist,
                                            sim_G11_real_A1.fdist,sim_G11_real_G2.fdist)),
                          mean.mhdist=mean(c(sim_A1_real_G2.mhdist,sim_A1_real_G11.mhdist,
                                             sim_G2_real_A1.mhdist,sim_G2_real_G11.mhdist,
                                             sim_G11_real_A1.mhdist,sim_G11_real_G2.mhdist)),
                          #
                          pair.mean.cor=mean(c(sim_A1_G2.cor$estimate, sim_A1_G11.cor$estimate, sim_G2_G11.cor$estimate)),
                          pair.mean.rho=mean(c(sim_A1_G2.rho$estimate, sim_A1_G11.rho$estimate, sim_G2_G11.rho$estimate)),
                          pair.mean.fdist=mean(c(sim_A1_G2.fdist, sim_A1_G11.fdist, sim_G2_G11.fdist)),
                          pair.mean.mhdist=mean(c(sim_A1_G2.mhdist, sim_A1_G11.mhdist, sim_G2_G11.mhdist))
        ))
      }) %>% bind_rows()  
    }, error = function(e){  
      cat("Error occurred in df_depth computation:", e$message, "\n")  
      NULL  # Return NULL if there's an error  
    })
    if (!is.null(df_depth_try)) {  
      break  # If df_depth computation is successful, break the loop  
    } else {  
      cat("Retrying due to error in df_depth computation...\n")  
    }
  }
  
  return(df_depth_try)
}) %>% bind_rows()

##output
out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/"
saveRDS(df_shuffle, file= paste0(out_path, "1000allDep_compare_shuffle_allLeafs_sim1.vs.real1",".Rds"))


## 3.2 compare with real transition rates (Figure 6B)
rm(list=ls())
path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig6_PLP_Lung/"
real.res <-readRDS(paste0(path, "compare_real.Rds"))
sim.res <- readRDS(paste0(path, "1000allDep_compare_shuffle_allLeafs_sim1.vs.real1.Rds"))

depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

#get p value base on rank
sim.real.res_p_rank <- lapply(seq(length(depth)), function(d){
  #d <- 5
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% 
    #select(depth, times, mean.cor,mean.fdist,mean.mhdist)
    select(depth, times, mean.cor,mean.rho,mean.fdist,mean.mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>%
    #select(depth, times, mean.cor,mean.fdist,mean.mhdist)
    select(depth, times, mean.cor,mean.rho,mean.fdist,mean.mhdist)
  sel.real.sim <- rbind(sel.op_shuffle, sel.op_real)
  col_name <- colnames(sel.real.sim)  
  #
  res_p <- mclapply(3:length(col_name), function(col){
    #col <- 3
    sel.col <- col_name[col]
    if (grepl("cor|rho", sel.col)) {
      sel.data <- sel.real.sim %>% select(c(1:2,sym(sel.col))) %>% .[order(.[,3], decreasing = T),]
    } else {
      sel.data <- sel.real.sim %>% select(c(1:2,sym(sel.col))) %>% .[order(.[,3], decreasing = F),]
    }
    p.value <- which(grepl("real", sel.data$times))/1000 #; 100-0.01
    return(data.frame(depth = sel.depth,
                      methods=sel.col,
                      p_value=p.value,
                      stringsAsFactors = F))
  }) %>% bind_rows()
  return(res_p)
}) %>% bind_rows() %>%
  filter(methods=="mean.fdist") %>% mutate(condition=paste0("<=",depth)) %>% select(-c(depth, methods))
colnames(sim.real.res_p_rank) <- c("p_rank", "condition")

#get p value base on test
sim.real.res_p.new <- lapply(seq(length(depth)), function(d){
  #d <- 5
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% 
    #select(depth, times, mean.cor,mean.fdist,mean.mhdist)
    select(depth, times, mean.cor,mean.rho,mean.fdist,mean.mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>%
    mutate(mean.cor=(A1_G2.cor+A1_G11.cor+G2_G11.cor)/3, sd.cor=sd(c(A1_G2.cor,A1_G11.cor,G2_G11.cor)),
           mean.rho=(A1_G2.rho+A1_G11.rho+G2_G11.rho)/3, sd.rho=sd(c(A1_G2.rho,A1_G11.rho,G2_G11.rho)),
           mean.fdist=(A1_G2.fdist+A1_G11.fdist+G2_G11.fdist)/3,
           sd.fdist=sd(c(A1_G2.fdist,A1_G11.fdist,G2_G11.fdist)), se.fdist=sd.fdist/sqrt(3),
           mean.mhdist=(A1_G2.mhdist+A1_G11.mhdist+G2_G11.mhdist)/3, sd.mhdist=sd(c(A1_G2.mhdist,A1_G11.mhdist,G2_G11.mhdist)))
  ##
  diff <- sel.op_real$mean.fdist-mean(sel.op_shuffle$mean.fdist)
  diff.mdist <- sel.op_real$mean.mhdist-mean(sel.op_shuffle$mean.mhdist)
  z_score <- (sel.op_real$mean.fdist - mean(sel.op_shuffle$mean.fdist)) / sd(sel.op_shuffle$mean.fdist)
  
  ##
  cor.p <- wilcox.test(sel.op_shuffle$mean.cor, mu=sel.op_real$mean.cor)$p.value %>% -log10(.)
  fdist.p <- wilcox.test(sel.op_shuffle$mean.fdist, mu=sel.op_real$mean.fdist)$p.value %>% -log10(.)
  mhdist.p <- wilcox.test(sel.op_shuffle$mean.mhdist, mu=sel.op_real$mean.mhdist)$p.value %>% -log10(.)
  
  return(data.frame(depth = sel.depth, Diff = diff, Diff.mdist = diff.mdist, z_score = z_score,
                    real.cor=sel.op_real$mean.cor, cor.Pvalue=cor.p,
                    real.fdist=sel.op_real$mean.fdist, real.fdist.sd=sel.op_real$sd.fdist, 
                    real.fdist.se=sel.op_real$se.fdist, fdist.Pvalue=fdist.p,
                    real.mhdist=sel.op_real$mean.mhdist, mhdist.Pvalue=mhdist.p,
                    stringsAsFactors = F))
}) %>% bind_rows() %>%
  mutate(condition=paste0("<=",depth)) %>% 
  #select(condition, Diff, z_score, fdist.Pvalue, real.fdist, real.fdist.sd, real.fdist.se)
  select(condition, z_score, fdist.Pvalue, real.fdist, real.fdist.sd, real.fdist.se)
##
sim.real.res_p.new <- full_join(sim.real.res_p.new, sim.real.res_p_rank, by="condition")


pdf("/mnt/data5/disk/yangwj/Result_plots/Figure6/PASTRI/Figure6B.new1.pdf",width=5,height=3.5)
ggplot(sim.real.res_p.new, aes(x=condition,y=real.fdist,color=condition))+geom_point(size=3)+#color="#1D91C0"
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Repeatability of robust\nAverage distance inferred from biological repeats")+
  geom_errorbar(aes(ymin=real.fdist-real.fdist.se, ymax=real.fdist+real.fdist.se),width=0)+
  theme_classic()+#scale_y_continuous(expand=c(0,0))+
  #coord_cartesian(ylim = c(min(sim.real.res_p.new$real.fdist - sim.real.res_p.new$real.fdist.se), 
  #                         max(sim.real.res_p.new$real.fdist + sim.real.res_p.new$real.fdist.se) * 1.15))+
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.6", "#91C893", "#D6D6D6")) +
  theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
        axis.text=element_text(size=10,color="black"),axis.line=element_line(color="black"), 
        panel.spacing = unit(1, "lines"), legend.position = "none")
dev.off()


#====== Figure 8D PhyloVelo ==========
values <- c(1.961395972196733, 1.9990496038232097, 2.035288657173599)
mean_val <- mean(values)
sd_val <- sd(values)
se_val <- sd_val/sqrt(3)
pdf("/mnt/data5/disk/yangwj/Result_plots/Figure6/PASTRI/FigureS8D.pdf",width=6,height=5)
ggplot(sim.real.res_p.new, aes(x=condition,y=real.fdist,color=condition))+geom_point(size=6)+#color="#1D91C0"
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Repeatability of robust\nAverage distance inferred from biological repeats")+
  geom_errorbar(aes(ymin=real.fdist-real.fdist.se, ymax=real.fdist+real.fdist.se),width=0)+
  theme_classic()+geom_hline(yintercept=mean_val, color="#E41A1C", linetype="dashed", linewidth=0.6) +
  geom_errorbar(aes(x=5, ymin=mean_val-se_val, ymax=mean_val+se_val), color="#E41A1C", alpha=0.3,linetype="dashed",linewidth=0.3,width=0) +
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.6", "#91C893", "#D6D6D6")) +
  theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
        axis.text=element_text(size=10,color="black"),axis.line=element_line(color="black"), 
        panel.spacing = unit(1, "lines"), legend.position = "none")
dev.off()

#====== Figure 6B significance test ==========
res.plot <- lapply(seq(length(depth)), function(d){
  #d <- 5
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(depth, mean.fdist)
  op_shuffle.1 <- sel.op_shuffle %>% melt(., id="depth")
  
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(depth, mean.fdist)
  
  #set fact label name
  facet_labels <- sim.real.res_p.new %>% filter(depth==sel.depth) %>% 
    mutate(label=paste0(condition, " fdist ", " P = ", p_rank)) 
  facet_labels <- setNames(facet_labels$label, "mean.fdist") 
  
  #plot
  p <- ggplot(op_shuffle.1, aes(x=value))+geom_histogram(binwidth=0.005,fill="#92C5DE")+
    geom_vline(data = sel.op_real, aes(xintercept = mean.fdist), linetype = "dashed", color="#D73027")+
    facet_wrap(~ variable, scales = "free", labeller = as_labeller(facet_labels), nrow=1)+
    theme_classic()+scale_y_continuous(expand=c(0,0))+ labs(x="",y="Count")+ #theme_bw()+
    theme(axis.title=element_text(size=10,color="black"), axis.text=element_text(size=10,color="black"),
          axis.line=element_line(color="black"), strip.text=element_text(size=8,color="black"),
          #strip.background=element_rect(fill="#E6E6E6"),strip.text.x=element_text(margin=unit(rep(8,4),"pt")),
          panel.spacing = unit(1, "lines"))
  return(p)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure6/PASTRI/FigureS8B.pdf",width=3,height=3)
res.plot[[5]]
dev.off()



