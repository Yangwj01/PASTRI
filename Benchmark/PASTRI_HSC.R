rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  #library(tidyverse)
  library(PASTRI)
  library(dplyr)
  library(parallel)
  library(ape)
  library(ggtree)
  library(patchwork)
  library(reshape2)
  library(tidyr)
  library(ggplot2)
  library(ggbreak)
  library(ggmagnify)
  library(RColorBrewer)
})

## ============ Run ===========================================================================
#== 1 get optimal transition
#calculate_lca_depths
tree_path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/"
B7.node.depth <- calculate_lca_depths(paste0(tree_path, "HSC_B7_B/HSC_B7.nwk"))
E12.node.depth <- calculate_lca_depths(paste0(tree_path, "HSC_E12_B/HSC_E12.nwk"))
F12.node.depth <- calculate_lca_depths(paste0(tree_path, "HSC_F12_B/HSC_F12.nwk"))
save(B7.node.depth, E12.node.depth, F12.node.depth,
     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/HSC.Node.depth.Rda")

#load some files
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/HSC.Node.depth.Rda")
cellInfo <- read.csv("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.all.cell.infos.csv")

#get optimal trans
Sample <- c("B7", "E12", "F12") 
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

df_mrca.mheight <- lapply(seq(length(Sample)), function(S){
  #S <- 3
  sel.Sample <- Sample[S]
  print(sel.Sample)
  sel.Onenode.depth <- paste0(sel.Sample,".node.depth") %>% get()
  sel.sample.name <- paste0("HSC_",sel.Sample)
  ##get cell info
  allelesInfo_path <- paste0("/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/",
                             sel.sample.name,"_B/",sel.sample.name,".AllelesInfo.csv")
  alleinfo <- read.csv(allelesInfo_path)
  sel.cell <- subset(cellInfo, cellInfo$orig.ident==sel.sample.name)
  
  sel.cell <- merge(alleinfo, sel.cell, by="BC") %>% dplyr::select(nodeLabel, BC, celltype, cellNum)
  colnames(sel.cell) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
  #sel.cell %>% pull(celltype) %>% table()
  sel.cell %>% filter(cellNum==1) %>% pull(celltype) %>% table()
  
  ##
  cell_type <- c("C4","C5","C6","C7","C8")
  Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
  Bound_matrix[upper.tri(Bound_matrix)] <- 0
  
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
  out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/"
  saveRDS(df_depth, file= paste0(out_path, sel.Sample,".op.trans",".Rds"))
  print("Successful")
  return(df_depth)
})

#get optimal trans bootstrap
depth <- 0.5
n <- 1000
df_mrca.mheight_bootstrap <- lapply(seq(length(Sample)), function(S){
  #S <- 1
  sel.Sample <- Sample[S]
  print(sel.Sample)
  sel.Onenode.depth <- paste0(sel.Sample,".node.depth") %>% get() 
  sel.sample.name <- paste0("HSC_",sel.Sample)
  ##get cell info
  allelesInfo_path <- paste0("/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/",
                             sel.sample.name,"_B/",sel.sample.name,".AllelesInfo.csv")
  alleinfo <- read.csv(allelesInfo_path)
  sel.cell <- subset(cellInfo, cellInfo$orig.ident==sel.sample.name)
  
  sel.cell <- merge(alleinfo, sel.cell, by="BC") %>% dplyr::select(nodeLabel, BC, celltype, cellNum)
  colnames(sel.cell) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
  #sel.cell %>% pull(celltype) %>% table()
  sel.cell %>% filter(cellNum==1) %>% pull(celltype) %>% table()
  
  #bootstrap get raw cell pair
  sel.onecell <- sel.cell %>% filter(cellNum==1) %>% filter(celltype %in% c("C4","C5","C6","C7","C8"))
  sel.Onenode.depth <- sel.Onenode.depth %>% filter(lca_normalized_height <= depth) %>%
    filter(node1 %in% sel.onecell$nodeLabel & node2 %in% sel.onecell$nodeLabel)
  
  ##
  cell_type <- c("C4","C5","C6","C7","C8")
  Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
  Bound_matrix[upper.tri(Bound_matrix)] <- 0
  
  ##Infer optimal constrained transition matrix 
  df_bootstrap <- lapply(seq(n), function(i){
    #i <- 1
    sel.depth <- depth
    bootstrap_node_depth <- sel.Onenode.depth[sample(nrow(sel.Onenode.depth), replace = TRUE), ]
    # Call get_optimal_transition_matrix() with the specified parameters
    mrca.mheight_optimal_results_d <- get_optimal_transition_matrix(
      node_pair_depth = bootstrap_node_depth,
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
  df_bootstrap_new <- df_bootstrap %>% group_by(type_combined_new) %>% 
    mutate(bootstrap_sd=sd(norm_optimal_T)) %>% select(-norm_optimal_T) %>% distinct()
  print("Successful")
  return(df_bootstrap_new)
}) %>% bind_rows()

saveRDS(df_mrca.mheight_bootstrap, 
        file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/op.trans.bootstrap.Rds")


#== 2 compare correlation between samples
#rm(list=ls())
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
B7.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/B7.op.trans.Rds")
E12.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/E12.op.trans.Rds")
F12.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/F12.op.trans.Rds")
cor.res <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  sel.depth.1 <- paste0("d_", sel.depth)
  print(sel.depth.1)
  ##
  irreversible.type <- c(#"C4_C4","C4_C5","C4_C6","C4_C7","C4_C8",
    "C5_C5","C5_C6","C5_C7","C5_C8",
    "C6_C6","C6_C7","C6_C8",
    "C7_C7","C7_C8",
    "C8_C8")
  sel.B7.op <- B7.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>%
    select(type_combined_new, start_cell, end_cell, norm_optimal_T)
  colnames(sel.B7.op) <- c("type.combined.new","start.cell","end.cell","B7.norm.optimal.T")
  sel.E12.op <- E12.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>% 
    select(type_combined_new, norm_optimal_T)
  colnames(sel.E12.op) <- c("type.combined.new","E12.norm.optimal.T")
  sel.F12.op <- F12.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>% 
    select(type_combined_new, norm_optimal_T)
  colnames(sel.F12.op) <- c("type.combined.new","F12.norm.optimal.T")
  all <- left_join(left_join(sel.B7.op, sel.E12.op, by="type.combined.new"), sel.F12.op, by="type.combined.new") %>%
    filter(start.cell!="C4" & end.cell!="C4")
  ##compare
  B7_E12.cor <- cor.test(all$B7.norm.optimal.T, all$E12.norm.optimal.T, method = "pearson")
  B7_E12.rho <- cor.test(all$B7.norm.optimal.T, all$E12.norm.optimal.T, method = "spearman")
  B7_E12.fdist <- dist(rbind(all$B7.norm.optimal.T, all$E12.norm.optimal.T), method = "euclidean") %>% as.vector()
  B7_E12.mhdist <- dist(rbind(all$B7.norm.optimal.T, all$E12.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  B7_F12.cor <- cor.test(all$B7.norm.optimal.T, all$F12.norm.optimal.T, method = "pearson")
  B7_F12.rho <- cor.test(all$B7.norm.optimal.T, all$F12.norm.optimal.T, method = "spearman")
  B7_F12.fdist <- dist(rbind(all$B7.norm.optimal.T, all$F12.norm.optimal.T), method = "euclidean") %>% as.vector()
  B7_F12.mhdist <- dist(rbind(all$B7.norm.optimal.T, all$F12.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  E12_F12.cor <- cor.test(all$E12.norm.optimal.T, all$F12.norm.optimal.T, method = "pearson")
  E12_F12.rho <- cor.test(all$E12.norm.optimal.T, all$F12.norm.optimal.T, method = "spearman")
  E12_F12.fdist <- dist(rbind(all$E12.norm.optimal.T, all$F12.norm.optimal.T), method = "euclidean") %>% as.vector()
  E12_F12.mhdist <- dist(rbind(all$E12.norm.optimal.T, all$F12.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  ##
  return(data.frame(depth=sel.depth,
                    times="real", 
                    stringsAsFactors =F,
                    B7_E12.cor=B7_E12.cor$estimate, B7_E12.cor.pval=B7_E12.cor$p.value,
                    B7_E12.rho=B7_E12.rho$estimate, B7_E12.rho.pval=B7_E12.rho$p.value,
                    B7_E12.fdist=B7_E12.fdist, B7_E12.mhdist=B7_E12.mhdist,
                    
                    B7_F12.cor=B7_F12.cor$estimate, B7_F12.cor.pval=B7_F12.cor$p.value,
                    B7_F12.rho=B7_F12.rho$estimate, B7_F12.rho.pval=B7_F12.rho$p.value,
                    B7_F12.fdist=B7_F12.fdist, B7_F12.mhdist=B7_F12.mhdist,
                    
                    E12_F12.cor=E12_F12.cor$estimate, E12_F12.cor.pval=E12_F12.cor$p.value,
                    E12_F12.rho=E12_F12.rho$estimate, E12_F12.rho.pval=E12_F12.rho$p.value,
                    E12_F12.fdist=E12_F12.fdist, E12_F12.mhdist=E12_F12.mhdist))
}) %>% bind_rows() %>%
  mutate(mean.cor=(B7_E12.cor+B7_F12.cor+E12_F12.cor)/3,
         mean.rho=(B7_E12.rho+B7_F12.rho+E12_F12.rho)/3,
         mean.fdist=(B7_E12.fdist+B7_F12.fdist+E12_F12.fdist)/3,
         mean.mhdist=(B7_E12.mhdist+B7_F12.mhdist+E12_F12.mhdist)/3)

saveRDS(cor.res, "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/compare_real.Rds")

## Figure 5B
p.cor.res <- lapply(seq(length(depth)), function(d){
  #d <- 4
  sel.depth <- depth[d]
  sel.depth.1 <- paste0("d_", sel.depth)
  print(sel.depth.1)
  ##
  irreversible.type <- c(#"C4_C4","C4_C5","C4_C6","C4_C7","C4_C8",
    "C5_C5","C5_C6","C5_C7","C5_C8",
    "C6_C6","C6_C7","C6_C8",
    "C7_C7","C7_C8",
    "C8_C8")
  sel.B7.op <- B7.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>%
    select(type_combined_new, start_cell, end_cell, norm_optimal_T)
  colnames(sel.B7.op) <- c("type.combined.new","start.cell","end.cell","B7.norm.optimal.T")
  sel.E12.op <- E12.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>% 
    select(type_combined_new, norm_optimal_T)
  colnames(sel.E12.op) <- c("type.combined.new","E12.norm.optimal.T")
  sel.F12.op <- F12.op %>% filter(depth==sel.depth.1) %>% filter(type_combined_new %in% irreversible.type) %>%
    select(type_combined_new, norm_optimal_T)
  colnames(sel.F12.op) <- c("type.combined.new","F12.norm.optimal.T")
  all <- left_join(left_join(sel.B7.op, sel.E12.op, by="type.combined.new"), sel.F12.op, by="type.combined.new")
  ##
  sel.cor.res <- cor.res %>% filter(depth==sel.depth)
  
  ##
  p.B7_E12 <- ggplot(all, aes(x = B7.norm.optimal.T, y = E12.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=3)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$B7_E12.cor, 2), " P=", round(sel.cor.res$B7_E12.cor.pval,10), 
                      #"\nrho = ", round(sel.cor.res$B7_E12.rho, 3), " P=", round(sel.cor.res$B7_E12.rho.pval,10),
                      "\neuclidean distance = ", round(sel.cor.res$B7_E12.fdist, 2)),
         x="PASTRI inferred transition rate of B7",y="PASTRI inferred transition rate of E12")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+xlim(0,1)+ylim(0,1)+
    scale_shape_manual(values=c(7,15,16,17))+scale_color_manual(values=c("#FF9116","#BBBFE0","#673A95","#B25A2C"))+
    theme(axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"),
          legend.position = "none")
  p.B7_F12 <- ggplot(all, aes(x = B7.norm.optimal.T, y = F12.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=3)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$B7_F12.cor, 2), " P=", round(sel.cor.res$B7_F12.cor.pval,10),  
                      #"\nrho = ", round(sel.cor.res$B7_F12.rho, 3), " P=", round(sel.cor.res$B7_F12.rho.pval,10),
                      "\neuclidean distance = ", round(sel.cor.res$B7_F12.fdist, 2)),
         x="PASTRI inferred transition rate of B7",y="PASTRI inferred transition rate of F12")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+xlim(0,1)+ylim(0,1)+
    scale_shape_manual(values=c(7,15,16,17))+scale_color_manual(values=c("#FF9116","#BBBFE0","#673A95","#B25A2C"))+
    theme(axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"),
          legend.position = "none")
  p.E12_F12 <- ggplot(all, aes(x = E12.norm.optimal.T, y = F12.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=3)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$E12_F12.cor, 2), " P=", round(sel.cor.res$E12_F12.cor.pval,10), 
                      #"\nrho = ", round(sel.cor.res$E12_F12.rho, 3), " P=", round(sel.cor.res$E12_F12.rho.pval,10),
                      "\neuclidean distance = ", round(sel.cor.res$E12_F12.fdist, 2)),
         x="PASTRI inferred transition rate of E12",y="PASTRI inferred transition rate of F12")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+xlim(0,1)+ylim(0,1)+
    scale_shape_manual(values=c(7,15,16,17))+scale_color_manual(values=c("#FF9116","#BBBFE0","#673A95","#B25A2C"))+
    theme(axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"),
          legend.text=element_text(size=10,color="black"))
  p <- p.B7_E12+p.B7_F12+p.E12_F12
  return(p)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/Figure5B.pdf",width=10,height=4)
p.cor.res[[4]]
dev.off()

#== 3 PASTRI accuracy
## 3.1 shuffle cellInfo and re-calculate transition rate by PASTRI
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
sel.sample <- "B7"
other.sample.1 <- "E12"
other.sample.2 <- "F12"
n <- 1000

df_shuffle <- lapply(seq(n), function(i){
  #i <- 1
  repeat{
    cat("times", i, "\n")
    sel.Onenode.depth <- paste0(sel.sample,".node.depth") %>% get()
    other1.Onenode.depth <- paste0(other.sample.1,".node.depth") %>% get()
    other2.Onenode.depth <- paste0(other.sample.2,".node.depth") %>% get()
    ##load sel sample celltype infos and random tree leaf
    sel.celltype <- subset(cellInfo, cellInfo$orig.ident==paste0("HSC_", sel.sample))
    other1.celltype <- subset(cellInfo, cellInfo$orig.ident==paste0("HSC_", other.sample.1))
    other2.celltype <- subset(cellInfo, cellInfo$orig.ident==paste0("HSC_", other.sample.2))
    #
    sel.allelesInfo <- read.csv(paste0("/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_",
                                       sel.sample,"_B/HSC_",sel.sample,".AllelesInfo.csv"))
    sel.raw.celltype <- merge(sel.allelesInfo, sel.celltype,by="BC") %>% dplyr::select(nodeLabel,BC,celltype,cellNum)
    colnames(sel.raw.celltype) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
    
    other1.allelesInfo <- read.csv(paste0("/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_",
                                          other.sample.1,"_B/HSC_",other.sample.1,".AllelesInfo.csv"))
    other1.raw.celltype <- merge(other1.allelesInfo, other1.celltype,by="BC") %>% dplyr::select(nodeLabel,BC,celltype,cellNum)
    colnames(other1.raw.celltype) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
    
    other2.allelesInfo <- read.csv(paste0("/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_",
                                          other.sample.2,"_B/HSC_",other.sample.2,".AllelesInfo.csv"))
    other2.raw.celltype <- merge(other2.allelesInfo, other2.celltype,by="BC") %>% dplyr::select(nodeLabel,BC,celltype,cellNum)
    colnames(other2.raw.celltype) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
    
    ##random tree leaf
    sel.sim.One.celltype <- sel.raw.celltype %>% filter(cellNum==1) %>% 
      #mutate(celltype=sample(celltype, size=nrow(.), replace=F))
      mutate(celltype=sample(sel.raw.celltype$celltype, size=nrow(.), replace=F))
    #table(sel.sim.One.celltype$celltype)
    other1.sim.One.celltype <- other1.raw.celltype %>% filter(cellNum==1) %>% 
      #mutate(celltype=sample(celltype, size=nrow(.), replace=F))
      mutate(celltype=sample(other1.raw.celltype$celltype, size=nrow(.), replace=F))
    other2.sim.One.celltype <- other2.raw.celltype %>% filter(cellNum==1) %>% 
      #mutate(celltype=sample(celltype, size=nrow(.), replace=F))
      mutate(celltype=sample(other2.raw.celltype$celltype, size=nrow(.), replace=F))
    #
    cell_type <- c("C4","C5","C6","C7","C8")
    Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
    Bound_matrix[upper.tri(Bound_matrix)] <- 0
    
    # Check if conditions are met
    df_depth_try <- tryCatch({  
      lapply(seq(length(depth)), function(d){
        #d <- 2
        sel.depth <- depth[d]
        print(sel.depth)
        ##get optimal transition matrix in shuffle trees by PASTRI
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
        B7.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/B7.op.trans.Rds")
        E12.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/E12.op.trans.Rds")
        F12.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/F12.op.trans.Rds")
        sel_real_B7 <- B7.op %>% filter(depth==paste0("d_",sel.depth)) %>% mutate(sample="real_B7")
        sel_real_E12 <- E12.op %>% filter(depth==paste0("d_",sel.depth)) %>% mutate(sample="real_E12")
        sel_real_F12 <- F12.op %>% filter(depth==paste0("d_",sel.depth)) %>% mutate(sample="real_F12")
        
        ##compare correlation between samples
        ##
        irreversible.type <- c(#"C4_C4","C4_C5","C4_C6","C4_C7","C4_C8",
          "C5_C5","C5_C6","C5_C7","C5_C8",
          "C6_C6","C6_C7","C6_C8",
          "C7_C7","C7_C8",
          "C8_C8")
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
        B7.op <- sel_real_B7 %>% filter(type_combined_new %in% irreversible.type) %>% select(type_combined_new, norm_optimal_T)
        colnames(B7.op) <- c("type.combined.new", "real.B7.norm.optimal.T")
        E12.op <- sel_real_E12 %>% filter(type_combined_new %in% irreversible.type) %>% select(type_combined_new, norm_optimal_T)
        colnames(E12.op) <- c("type.combined.new", "real.E12.norm.optimal.T")
        F12.op <- sel_real_F12 %>% filter(type_combined_new %in% irreversible.type) %>% select(type_combined_new, norm_optimal_T)
        colnames(F12.op) <- c("type.combined.new", "real.F12.norm.optimal.T")
        
        all <- full_join(full_join(full_join(sel.op, other_1.op, by="type.combined.new"),
                                   full_join(other_2.op, B7.op, by="type.combined.new"), by="type.combined.new"),
                         full_join(E12.op, F12.op, by="type.combined.new"), by="type.combined.new")
        ##compare
        sim_B7_real_E12.cor <- cor.test(all$B7.norm.optimal.T, all$real.E12.norm.optimal.T, method = "pearson")
        sim_B7_real_E12.rho <- cor.test(all$B7.norm.optimal.T, all$real.E12.norm.optimal.T, method = "spearman")
        sim_B7_real_E12.fdist <- dist(rbind(all$B7.norm.optimal.T, all$real.E12.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_B7_real_E12.mhdist <- dist(rbind(all$B7.norm.optimal.T, all$real.E12.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_B7_real_F12.cor <- cor.test(all$B7.norm.optimal.T, all$real.F12.norm.optimal.T, method = "pearson")
        sim_B7_real_F12.rho <- cor.test(all$B7.norm.optimal.T, all$real.F12.norm.optimal.T, method = "spearman")
        sim_B7_real_F12.fdist <- dist(rbind(all$B7.norm.optimal.T, all$real.F12.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_B7_real_F12.mhdist <- dist(rbind(all$B7.norm.optimal.T, all$real.F12.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_E12_real_B7.cor <- cor.test(all$E12.norm.optimal.T, all$real.B7.norm.optimal.T, method = "pearson")
        sim_E12_real_B7.rho <- cor.test(all$E12.norm.optimal.T, all$real.B7.norm.optimal.T, method = "spearman")
        sim_E12_real_B7.fdist <- dist(rbind(all$E12.norm.optimal.T, all$real.B7.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_E12_real_B7.mhdist <- dist(rbind(all$E12.norm.optimal.T, all$real.B7.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_E12_real_F12.cor <- cor.test(all$E12.norm.optimal.T, all$real.F12.norm.optimal.T, method = "pearson")
        sim_E12_real_F12.rho <- cor.test(all$E12.norm.optimal.T, all$real.F12.norm.optimal.T, method = "spearman")
        sim_E12_real_F12.fdist <- dist(rbind(all$E12.norm.optimal.T, all$real.F12.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_E12_real_F12.mhdist <- dist(rbind(all$E12.norm.optimal.T, all$real.F12.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_F12_real_B7.cor <- cor.test(all$F12.norm.optimal.T, all$real.B7.norm.optimal.T, method = "pearson")
        sim_F12_real_B7.rho <- cor.test(all$F12.norm.optimal.T, all$real.B7.norm.optimal.T, method = "spearman")
        sim_F12_real_B7.fdist <- dist(rbind(all$F12.norm.optimal.T, all$real.B7.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_F12_real_B7.mhdist <- dist(rbind(all$F12.norm.optimal.T, all$real.B7.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_F12_real_E12.cor <- cor.test(all$F12.norm.optimal.T, all$real.E12.norm.optimal.T, method = "pearson")
        sim_F12_real_E12.rho <- cor.test(all$F12.norm.optimal.T, all$real.E12.norm.optimal.T, method = "spearman")
        sim_F12_real_E12.fdist <- dist(rbind(all$F12.norm.optimal.T, all$real.E12.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_F12_real_E12.mhdist <- dist(rbind(all$F12.norm.optimal.T, all$real.E12.norm.optimal.T), method = "manhattan") %>% as.vector()
        ##
        return(data.frame(depth=sel.depth,
                          times=i, 
                          stringsAsFactors =F,
                          sim_B7_real_E12.cor=sim_B7_real_E12.cor$estimate, sim_B7_real_E12.cor.pval=sim_B7_real_E12.cor$p.value,
                          sim_B7_real_E12.rho=sim_B7_real_E12.rho$estimate, sim_B7_real_E12.rho.pval=sim_B7_real_E12.rho$p.value,
                          sim_B7_real_E12.fdist=sim_B7_real_E12.fdist, sim_B7_real_E12.mhdist=sim_B7_real_E12.mhdist,
                          
                          sim_B7_real_F12.cor=sim_B7_real_F12.cor$estimate, sim_B7_real_F12.cor.pval=sim_B7_real_F12.cor$p.value,
                          sim_B7_real_F12.rho=sim_B7_real_F12.rho$estimate, sim_B7_real_F12.rho.pval=sim_B7_real_F12.rho$p.value,
                          sim_B7_real_F12.fdist=sim_B7_real_F12.fdist, sim_B7_real_F12.mhdist=sim_B7_real_F12.mhdist,
                          
                          sim_E12_real_B7.cor=sim_E12_real_B7.cor$estimate, sim_E12_real_B7.cor.pval=sim_E12_real_B7.cor$p.value,
                          sim_E12_real_B7.rho=sim_E12_real_B7.rho$estimate, sim_E12_real_B7.rho.pval=sim_E12_real_B7.rho$p.value,
                          sim_E12_real_B7.fdist=sim_E12_real_B7.fdist, sim_E12_real_B7.mhdist=sim_E12_real_B7.mhdist,
                          
                          sim_E12_real_F12.cor=sim_E12_real_F12.cor$estimate, sim_E12_real_F12.cor.pval=sim_E12_real_F12.cor$p.value,
                          sim_E12_real_F12.rho=sim_E12_real_F12.rho$estimate, sim_E12_real_F12.rho.pval=sim_E12_real_F12.rho$p.value,
                          sim_E12_real_F12.fdist=sim_E12_real_F12.fdist, sim_E12_real_F12.mhdist=sim_E12_real_F12.mhdist,
                          
                          sim_F12_real_B7.cor=sim_F12_real_B7.cor$estimate, sim_F12_real_B7.cor.pval=sim_F12_real_B7.cor$p.value,
                          sim_F12_real_B7.rho=sim_F12_real_B7.rho$estimate, sim_F12_real_B7.rho.pval=sim_F12_real_B7.rho$p.value,
                          sim_F12_real_B7.fdist=sim_F12_real_B7.fdist, sim_F12_real_B7.mhdist=sim_F12_real_B7.mhdist,
                          
                          sim_F12_real_E12.cor=sim_F12_real_E12.cor$estimate, sim_F12_real_E12.cor.pval=sim_F12_real_E12.cor$p.value,
                          sim_F12_real_E12.rho=sim_F12_real_E12.rho$estimate, sim_F12_real_E12.rho.pval=sim_F12_real_E12.rho$p.value,
                          sim_F12_real_E12.fdist=sim_F12_real_E12.fdist, sim_F12_real_E12.mhdist=sim_F12_real_E12.mhdist,
                          
                          mean.cor=mean(c(sim_B7_real_E12.cor$estimate, sim_B7_real_F12.cor$estimate,
                                          sim_E12_real_B7.cor$estimate, sim_E12_real_F12.cor$estimate,
                                          sim_F12_real_B7.cor$estimate, sim_F12_real_E12.cor$estimate)),
                          mean.rho=mean(c(sim_B7_real_E12.rho$estimate, sim_B7_real_F12.rho$estimate,
                                          sim_E12_real_B7.rho$estimate, sim_E12_real_F12.rho$estimate,
                                          sim_F12_real_B7.rho$estimate, sim_F12_real_E12.rho$estimate)),
                          mean.fdist=mean(c(sim_B7_real_E12.fdist, sim_B7_real_F12.fdist,
                                            sim_E12_real_B7.fdist, sim_E12_real_F12.fdist,
                                            sim_F12_real_B7.fdist, sim_F12_real_E12.fdist)),
                          mean.mhdist=mean(c(sim_B7_real_E12.mhdist, sim_B7_real_F12.mhdist,
                                             sim_E12_real_B7.mhdist, sim_E12_real_F12.mhdist,
                                             sim_F12_real_B7.mhdist, sim_F12_real_E12.mhdist))
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
out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/"
saveRDS(df_shuffle, file= paste0(out_path, "1000allDep_compare_shuffle_allLeafs_sim1.vs.real1",".Rds"))


## 3.2 compare with real transition rates (Figure 5A)
rm(list=ls())
real.res <-readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/compare_real.Rds")
sim.res <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/1000allDep_compare_shuffle_allLeafs_sim1.vs.real1.Rds")

depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

sim.real.res_p.new <- lapply(seq(length(depth)), function(d){
  #d <- 5
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>% 
    select(depth, times, mean.cor, mean.fdist,mean.mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>%
    mutate(mean.cor=(B7_E12.cor+B7_F12.cor+E12_F12.cor)/3, sd.cor=sd(c(B7_E12.cor,B7_F12.cor,E12_F12.cor)),
           mean.fdist=(B7_E12.fdist+B7_F12.fdist+E12_F12.fdist)/3, 
           sd.fdist=sd(c(B7_E12.fdist,B7_F12.fdist,E12_F12.fdist)), se.fdist=sd.fdist/sqrt(3),
           mean.mhdist=(B7_E12.mhdist+B7_F12.mhdist+E12_F12.mhdist)/3,
           sd.mhdist=sd(c(B7_E12.mhdist, B7_F12.mhdist, E12_F12.mhdist)))
  ##
  diff <- sel.op_real$mean.fdist-mean(sel.op_shuffle$mean.fdist)
  z_score <- (sel.op_real$mean.fdist - mean(sel.op_shuffle$mean.fdist)) / sd(sel.op_shuffle$mean.fdist)
  
  ##
  cor.p <- wilcox.test(sel.op_shuffle$mean.cor, mu=sel.op_real$mean.cor)$p.value %>% -log10(.)
  fdist.p <- wilcox.test(sel.op_shuffle$mean.fdist, mu=sel.op_real$mean.fdist)$p.value %>% -log10(.)
  mhdist.p <- wilcox.test(sel.op_shuffle$mean.mhdist, mu=sel.op_real$mean.mhdist)$p.value %>% -log10(.)
  
  return(data.frame(depth = sel.depth, Diff = diff, z_score = z_score,
                    real.cor=sel.op_real$mean.cor, cor.Pvalue=cor.p,
                    real.fdist=sel.op_real$mean.fdist, real.fdist.sd=sel.op_real$sd.fdist, 
                    real.fdist.se=sel.op_real$se.fdist, fdist.Pvalue=fdist.p,
                    real.mhdist=sel.op_real$mean.mhdist, mhdist.Pvalue=mhdist.p,
                    stringsAsFactors = F))
}) %>% bind_rows() %>%
  mutate(condition=paste0("<=",depth)) %>% 
  #select(condition, Diff, z_score, fdist.Pvalue, real.fdist, real.fdist.sd, real.fdist.se)
  select(condition, z_score, fdist.Pvalue, real.fdist, real.fdist.sd, real.fdist.se)

#====== Figure 5A ==========
#get p value base on rank
sim.real.res_p_rank <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>% 
    select(depth, times, mean.cor,mean.rho,mean.fdist,mean.mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>%
    mutate(mean.cor=(B7_E12.cor+B7_F12.cor+E12_F12.cor)/3,
           mean.rho=(B7_E12.rho+B7_F12.rho+E12_F12.rho)/3,
           mean.fdist=(B7_E12.fdist+B7_F12.fdist+E12_F12.fdist)/3,
           mean.mhdist=(B7_E12.mhdist+B7_F12.mhdist+E12_F12.mhdist)/3) %>% 
    select(depth, times, mean.cor, mean.rho, mean.fdist,mean.mhdist)
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
sim.real.res_p_rank

##
sim.real.res_p.new <- full_join(sim.real.res_p.new, sim.real.res_p_rank, by="condition")

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/Figure5A.new1.pdf",width=4.8,height=3.5)
ggplot(sim.real.res_p.new, aes(x = condition, y = real.fdist, color=condition)) +geom_point(size=3)+
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Repeatability of robust\nAverage distance inferred from biological repeats")+
  geom_errorbar(aes(ymin=real.fdist-real.fdist.se, ymax=real.fdist+real.fdist.se),width=0)+
  theme_classic()+#scale_y_continuous(expand=c(0,0))+
  #coord_cartesian(ylim = c(min(sim.real.res_p.new$real.fdist - sim.real.res_p.new$real.fdist.se), 
  #                         max(sim.real.res_p.new$real.fdist + sim.real.res_p.new$real.fdist.se) * 1.15))+
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.5", "#91C893", "#D6D6D6"))+
  theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
        axis.text=element_text(size=10,color="black"),axis.line=element_line(color="black"), 
        panel.spacing = unit(1, "lines"), legend.position = "none")
dev.off()


#====== Figure S8A ==========
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

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/FigureS8A.pdf",width=3,height=3)
res.plot[[4]]
dev.off()

#====== Figure S8B PhyloVelo ==========
values <- c(0.71628862782321, 1.5429902547657461, 1.3169737824652772)
mean_val <- mean(values)
sd_val <- sd(values)
se_val <- sd_val/sqrt(3)
pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/FigureS8B.pdf",width=6,height=5)
ggplot(sim.real.res_p.new, aes(x = condition, y = real.fdist, color=condition)) +geom_point(size=6)+
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Repeatability of robust\nAverage distance inferred from biological repeats")+
  geom_errorbar(aes(ymin=real.fdist-real.fdist.se, ymax=real.fdist+real.fdist.se),width=0)+
  theme_classic()+geom_hline(yintercept=mean_val, color="#E41A1C", linetype="dashed", linewidth=0.6) +
  geom_errorbar(aes(x=5, ymin=mean_val-se_val, ymax=mean_val+se_val), color="#E41A1C", alpha=0.3,linetype="dashed",,width=0) +
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.5", "#91C893", "#D6D6D6"))+
  theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
        axis.text=element_text(size=10,color="black"),axis.line=element_line(color="black"), 
        panel.spacing = unit(1, "lines"), legend.position = "none")
dev.off()


#====== Fig 5F ======
#compare 0.5<d<=0.6 by PASTRI
rm(list=ls())
B7.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/B7.op.trans.Rds") %>%
  filter(depth=="d_0.5") %>% select(-c(depth,sample))
colnames(B7.op) <- c("type.combined.new","start.cell","end.cell","B7")
E12.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/E12.op.trans.Rds") %>%
  filter(depth=="d_0.5") %>% select(type_combined_new,norm_optimal_T)
colnames(E12.op) <- c("type.combined.new","E12")
F12.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/F12.op.trans.Rds") %>%
  filter(depth=="d_0.5") %>% select(type_combined_new,norm_optimal_T)
colnames(F12.op) <- c("type.combined.new","F12")
all <- full_join(full_join(B7.op, E12.op, by="type.combined.new"), F12.op, by="type.combined.new") %>%
  mutate(mean_value=apply(.[,c("B7", "E12", "F12")], 1, mean),
         sd_value=apply(.[,c("B7", "E12", "F12")], 1, sd)) %>% #filter(mean_value>0)
  filter(type.combined.new %in% c("C4_C4","C4_C5","C4_C6","C4_C7","C4_C8",
                                  "C5_C5","C5_C6","C5_C7","C5_C8",
                                  "C6_C6","C6_C7","C6_C8",
                                  "C7_C7","C7_C8","C8_C8")) 

all.new <- all %>% filter(start.cell!="C4") %>% mutate(type=ifelse(start.cell=="C5", "Early", "Late"))
all.range <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/range.all.op.trans.Rds") %>%
  filter(start.cell!="C4") %>% mutate(type=ifelse(start.cell=="C5", "Early", "Late"))
colnames(all.range) <- c("type.combined.new","start.cell","end.cell","range.B7","range.E12","range.F12","range.mean_value","range.sd_value","type")
d1_d2 <- full_join(all.new, all.range, by=c("type.combined.new","start.cell","end.cell","type")) %>%
  mutate(diff_B7=range.B7-B7, diff_E12=range.E12-E12, diff_F12=range.F12-F12, 
         diff_mean=range.mean_value-mean_value, diff_sd=range.sd_value-sd_value,
         diff_B7_ratio=diff_B7/B7, diff_E12_ratio=diff_E12/E12, diff_F12_ratio=diff_F12/F12,
         diff_mean_ratio=diff_mean/mean_value, diff_sd_ratio=diff_sd/sd_value) %>% select(c(1,9,15:24))


d1_d2.new <- d1_d2[,c(1:2, 6:7)] %>% filter(!type.combined.new %in% c("C7_C7","C7_C8","C8_C8"))
d1_d2.new$type.combined.new <- factor(d1_d2.new$type.combined.new, levels = c("C5_C5","C5_C6","C6_C6",
                                                                              "C5_C7","C6_C7","C5_C8","C6_C8"))
#d1_d2.new$type <- factor(d1_d2.new$type, levels = c("Early", "Late"))

p1 <- ggplot(d1_d2.new, aes(x=type.combined.new, y=diff_mean))+ 
  geom_bar(stat="identity", position="identity", width=0.6, fill="#4C8CB8", alpha=0.36)+
  geom_hline(yintercept=0, color="gray60", linetype = "dashed")+
  labs(x = "", y = "Relative deviation of\n inferred transition rate\n(0.5 < d ≤ 0.6 v.s. d ≤ 0.5)")+ theme_classic()+
  theme(axis.text.x=element_text(color="black", size=10, angle=30,hjust = 1,vjust = 1), 
        axis.text.y=element_text(color="black",size=10), axis.title.y=element_text(color="black", size=10),
        axis.line=element_line(color="black"))
p2 <- ggplot(d1_d2.new, aes(x=type.combined.new, y=diff_sd))+ 
  geom_bar(stat="identity", position="identity", width=0.6, fill="#4C8CB8", alpha=0.36)+
  geom_hline(yintercept=0, color="gray60", linetype = "dashed")+
  labs(x = "", y = "Decrease in reproducibility\nof transition rate inference\n(0.5 < d ≤ 0.6 v.s. d ≤ 0.5)")+ theme_classic()+
  theme(axis.text.x=element_text(color="black", size=10, angle=30,hjust = 1,vjust = 1), 
        axis.text.y=element_text(color="black",size=10), axis.title.y=element_text(color="black", size=10),
        axis.line=element_line(color="black"))

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/Figure5F.pdf",width=3.5,height=4.5)
p1/p2
dev.off()


#====================================================================================================
##Calculate the transcriptional differences among different cell types
rm(list=ls())
library(Seurat)
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")
unique(sce$orig.ident)

sel.depth = 0.5
Sample <- c("B7", "E12", "F12")
celltype_list=c("C4","C5","C6","C7","C8")
exp.ratio <- 0.5

Edist.allcell <- lapply(seq(length(Sample)), function(s){
  #s <- 1
  sel.Sample <- Sample[s]
  print(sel.Sample)
  sel.sce <- subset(sce, subset = orig.ident == paste0("10X_HSC_", sel.Sample))
  sel.data <- sel.sce@meta.data %>% select(celltype) %>% filter(celltype %in% celltype_list)
  print(table(sel.data$celltype))
  
  # 1 enumerate all cell pairs
  cells <- rownames(sel.data)
  cell_pairs_combinations <- combn(unique(cells), 2)
  cell_pairs_df <- data.frame(
    node1 = cell_pairs_combinations[1, ],
    node2 = cell_pairs_combinations[2, ],
    stringsAsFactors = FALSE
  )
  cell_pairs_df$node1_type <- sel.data$celltype[match(cell_pairs_df$node1, rownames(sel.data))]
  cell_pairs_df$node2_type <- sel.data$celltype[match(cell_pairs_df$node2, rownames(sel.data))]
  
  # 2.1 Construct standardized cell type combinations to avoid duplicate counts
  unique_type_combn <- t(combn(unique(celltype_list), 2)) %>% as.data.frame(stringsAsFactors = FALSE)
  colnames(unique_type_combn) <- c("type_1", "type_2")
  unique_type_combn <- bind_rows(
    data.frame(type_1 = celltype_list, type_2 = celltype_list, stringsAsFactors = FALSE),
    unique_type_combn
  )
  unique_type_combn$ref_comb <- paste(unique_type_combn$type_1, unique_type_combn$type_2, sep = "_")
  unique_type_combn$reverse_comb <- paste(unique_type_combn$type_2, unique_type_combn$type_1, sep = "_")
  
  combn_standard_df <- data.frame(
    ref_comb = c(unique_type_combn$ref_comb, unique_type_combn$ref_comb),
    present_comb = c(unique_type_combn$ref_comb, unique_type_combn$reverse_comb),
    stringsAsFactors = FALSE
  )
  rm(unique_type_combn)
  
  # 2.2 Standardize the cell type combination in node_pair_depth
  cell_pairs_df$type_combined <- paste(cell_pairs_df$node1_type, cell_pairs_df$node2_type, sep = "_")
  cell_pairs_df$type_combined <- combn_standard_df$ref_comb[match(cell_pairs_df$type_combined, combn_standard_df$present_comb)]
  rm(combn_standard_df)
  
  # 3 get exp matrix, filter more than 50% cell exp genes=
  exp.matrix <- sel.sce@assays$RNA@data %>% as.data.frame()
  exp.matrix$freq <- apply(exp.matrix, 1, function(x){sum(x>0)/ncol(exp.matrix)})
  exp.matrix <- exp.matrix[exp.matrix$freq > exp.ratio,] %>% select(-ncol(exp.matrix))
  print(nrow(exp.matrix))
  
  ExpDist <- mclapply(seq(nrow(cell_pairs_df)), function(j){
    #j <- 1
    node1 <- cell_pairs_df$node1[j]
    node2 <- cell_pairs_df$node2[j]
    # get exp 
    node1.exp <- exp.matrix[, node1] %>% as.numeric()
    node2.exp <- exp.matrix[, node2] %>% as.numeric()
    #expDist <- dist(rbind(node1.exp, node2.exp), method = "euclidean")[1]
    #expDist2 <- 1-cor(node1.exp, node2.exp,method="pearson")
    expDist <- amap::Dist(rbind(node1.exp,node2.exp),method="euclidean")[1]
    expDist2 <- amap::Dist(rbind(node1.exp,node2.exp),method="pearson")[1]
    return(data.frame(node1=node1,
                      node2=node2,
                      node1_celltype=cell_pairs_df$node1_type[j],
                      node2_celltype=cell_pairs_df$node2_type[j],
                      type_combined=cell_pairs_df$type_combined[j],
                      expDist.eu=expDist,
                      expDist.cor=expDist2,
                      stringsAsFactors = F))
  }, mc.cores = 60) %>% bind_rows()
  #unique(ExpDist$type_combined)
  ExpDist.new <- ExpDist %>% ungroup() %>% group_by(type_combined) %>% 
    mutate(celltype.ExpDist.mean.eu=mean(expDist.eu), celltype.ExpDist.mean.cor=mean(expDist.cor)) %>%
    ungroup() %>% select(type_combined, celltype.ExpDist.mean.eu, celltype.ExpDist.mean.cor) %>% unique()
  
  # 4 get transition matrix at D <= 0.5
  sel.op <- readRDS(paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/",
                           sel.Sample, ".op.trans.Rds")) %>% filter(depth==paste0("d_", sel.depth))
  # 5 compare ExpDist with transition rate
  Edist.op <- merge(ExpDist.new, sel.op, by.x="type_combined", by.y="type_combined_new") %>% 
    select(type_combined, start_cell, end_cell, norm_optimal_T, celltype.ExpDist.mean.eu, celltype.ExpDist.mean.cor)
  colnames(Edist.op) <- c("type_combined", "start_cell", "end_cell", paste0(sel.Sample, ".T"), paste0(sel.Sample, ".eu"), paste0(sel.Sample, ".cor"))
  return(Edist.op)
  #cor.eu <- cor.test(Edist.op$norm_optimal_T, Edist.op$celltype.ExpDist.mean.eu, method = "pearson")
  #cor.cor <- cor.test(Edist.op$norm_optimal_T, Edist.op$celltype.ExpDist.mean.cor, method = "pearson")
  
  #res <- data.frame(sample = sel.Sample,
  #                  cor.Expeu.vs.L = cor.eu$estimate, cor.Expeu.vs.L.P = cor.eu$p.value,
  #                  cor.ExpPearson.vs.L = cor.cor$estimate, cor.ExpPearson.vs.L.P = cor.cor$p.value,
  #                  stringsAsFactors = F)
  #return(res)
})

Edist.allcell.new <- full_join(full_join(Edist.allcell[[1]], Edist.allcell[[2]], by=c("type_combined", "start_cell", "end_cell")),
                               Edist.allcell[[3]], by=c("type_combined", "start_cell", "end_cell")) %>%
  filter(start_cell != "C4") %>% #filter(!type_combined %in% c("C8_C8")) %>%
  mutate(mean.T=(B7.T+E12.T+F12.T)/3, 
         mean.eu=(B7.eu+E12.eu+F12.eu)/3, mean.cor=(B7.cor+E12.cor+F12.cor)/3)
saveRDS(Edist.allcell.new, file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/ExpDiff_vs_TransRate_expRatio0.5.Rds")

#====
rm(list=ls())
path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig5_HSC/"
Edist.allcell.new <- readRDS(paste0(path, "ExpDiff_vs_TransRate_expRatio0.5.Rds"))

Edist.allcell.new <- Edist.allcell.new[,c(1:3,13:15)]
real.res <- data.frame(cor.eu = as.numeric(cor.test(Edist.allcell.new$mean.T, Edist.allcell.new$mean.eu, method="spearman")$estimate),
                       cor.eu.P = cor.test(Edist.allcell.new$mean.T, Edist.allcell.new$mean.eu, method="spearman")$p.value,
                       cor.pearson = as.numeric(cor.test(Edist.allcell.new$mean.T, Edist.allcell.new$mean.cor, method="spearman")$estimate),
                       cor.pearson.P = cor.test(Edist.allcell.new$mean.T, Edist.allcell.new$mean.cor, method="spearman")$p.value,
                       times = "real",
                       stringsAsFactors = F)
real.res
##suffle 4!*3!*2!*1!=288
library(purrr)
generate_all_permutations_simple <- function(df) {
  #df <- Edist.allcell.new
  result <- df %>%
    group_by(start_cell) %>%
    group_modify(~ {
      group_values <- .x$mean.T
      n <- length(group_values)
      
      if(n == 1) {
        permutations <- list(group_values)
      } else {
        permutations <- combinat::permn(group_values)
      }
      
      map_dfr(seq_along(permutations), function(idx) {
        .x %>%
          mutate(
            T_shuffled = permutations[[idx]],
            permutation_within_group = idx,
            row_id = row_number()
          )
      })
    }) %>%
    ungroup()
  
  # combine all the permutations across groups
  # Obtain the permutation IDs for each group
  group_perm_ids <- result %>%
    group_by(start_cell) %>%
    summarise(
      perm_ids = list(unique(permutation_within_group)),
      .groups = 'drop'
    )
  
  # Generate all possible combinations
  args_list <- map2(group_perm_ids$start_cell,
                    group_perm_ids$perm_ids,
                    ~ setNames(list(.y), .x)) %>% unlist(recursive = FALSE)
  all_combinations <- expand.grid(args_list) %>% mutate(combination_id = row_number())
  
  # Reconstruct the data based on the combination ID
  final_data <- map_dfr(1:nrow(all_combinations), function(i) {
    #i <- 1
    comb <- all_combinations[i, ]
    
    # Select the corresponding arrangement from each group
    bind_rows(
      result %>% filter(start_cell == "C5", permutation_within_group == comb$C5),
      result %>% filter(start_cell == "C6", permutation_within_group == comb$C6),
      result %>% filter(start_cell == "C7", permutation_within_group == comb$C7),
      result %>% filter(start_cell == "C8", permutation_within_group == comb$C8)
    ) %>%
      mutate(combination_id = i)
  })
  
  return(final_data)
}

# Generate all permutations and combinations
all_permutations_df <- generate_all_permutations_simple(Edist.allcell.new)

times <- max(all_permutations_df$combination_id)
sim.res <- lapply(seq(times), function(t){
  #t <- 1
  sel.data <- all_permutations_df %>% filter(combination_id==t)
  sim.res1 <- data.frame(cor.eu = as.numeric(cor.test(sel.data$T_shuffled, sel.data$mean.eu, method="spearman")$estimate),
                         cor.eu.P = cor.test(sel.data$T_shuffled, sel.data$mean.eu, method="spearman")$p.value,
                         cor.pearson = as.numeric(cor.test(sel.data$T_shuffled, sel.data$mean.cor, method="spearman")$estimate),
                         cor.pearson.P = cor.test(sel.data$T_shuffled, sel.data$mean.cor, method="spearman")$p.value,
                         times = t,
                         stringsAsFactors = F)
  return(sim.res1)
}) %>% bind_rows()

col_name <- colnames(sim.res[,c(1,3)]) 
res_p <- lapply(seq(length(col_name)), function(col){
  #col <- 1
  sel.col <- col_name[col]
  sel.data <- rbind(real.res, sim.res) %>% select(c(sym(sel.col), 5)) %>% .[order(.[,1], decreasing = F),]
  p.value <- which(grepl("real", sel.data$times))/287 #; 100-0.01
  return(data.frame(method=sel.col,
                    p_value=p.value,
                    stringsAsFactors = F))
}) %>% bind_rows()
res_p

max_density <- max(ggplot_build(ggplot(sim.res, aes(x = cor.pearson)) +  geom_histogram(binwidth=0.03,fill="#92C5DE"))$data[[1]]$count)


pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/Figure5D.pdf",width=3.2,height=3)
ggplot(sim.res, aes(x=cor.pearson))+ geom_histogram(binwidth=0.03,fill="#92C5DE", alpha=0.6)+
  #geom_density(fill="#92C5DE", alpha=0.5)+
  geom_segment(data=real.res, 
               aes(x=cor.pearson, xend=cor.pearson ,y=max_density*0.6, yend=0, color="red"),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"))+
  labs(x="Transcriptomic difference vs transition rate",y="Count")+
  theme_classic()+scale_y_continuous(limits=c(0, max_density*1.1), expand=c(0,0))+
  theme(axis.title=element_text(size=10,color="black"), axis.text=element_text(size=10,color="black"),
        legend.position = "none",axis.line=element_line(color="black"))
dev.off()


#=======================================================================================================
#transition fast, the fewer the differentially expressed genes
## ====== import some librarys ======
rm(list=ls())
library(Seurat)
library(dplyr)

## ====== source and define functions ======
find_filtered_markers <- function(sce_obj, group1, group2, filter_genes, min_pct = 0.25, padj_thresh = 0.05) {
  markers <- FindMarkers(sce_obj, ident.1 = group1, ident.2 = group2, min.pct = min_pct)
  markers_filtered <- markers %>% filter(p_val_adj < padj_thresh) %>% filter(rownames(.) %in% filter_genes)
  
  markers_filtered$contrast <- paste0(group1, "_", group2)
  markers_filtered$gene <- rownames(markers_filtered)
  
  return(markers_filtered)
}

## ====== Run ======
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")
sce <- subset(sce, subset = orig.ident %in% c("10X_HSC_B7", "10X_HSC_E12", "10X_HSC_F12"))
Idents(sce) <- "celltype"

# get expression matrix
exp.ratio = 0.5
exp.matrix <- sce@assays$RNA@data %>% as.data.frame()
exp.matrix$freq <- apply(exp.matrix, 1, function(x){sum(x>0)/ncol(exp.matrix)})
exp.matrix <- exp.matrix[exp.matrix$freq > exp.ratio,] %>% select(-ncol(exp.matrix))
filter_gene <- rownames(exp.matrix)


#compare
contrasts <- list(
  # more than mean rate 0.4
  c("C5", "C7"), c("C6", "C7"),
  # less than mean rate 0.4
  c("C5", "C6"), c("C5", "C8"), c("C6", "C8"), c("C7", "C8"))

res <- lapply(seq(length(contrasts)), function(c){
  #c <- 1
  sel.contrast <- contrasts[[c]]
  cat("contrast ", sel.contrast)
  marker <- find_filtered_markers(sce_obj = sce,
                                  group1=sel.contrast[1],
                                  group2=sel.contrast[2],
                                  filter_genes = filter_gene,
                                  min_pct = 0.25,
                                  padj_thresh = 0.05)
  return(marker)
}) %>% bind_rows()

res1 <- res %>% group_by(contrast) %>% mutate(DE_num=length(gene), No_DE_num=length(setdiff(filter_gene, gene))) %>%
  ungroup() %>% select(contrast,DE_num, No_DE_num) %>% distinct() %>% as.data.frame() %>%
  mutate(class = case_when(contrast %in% c("C5_C7", "C6_C7") ~ "more",
                           contrast %in% c("C5_C6", "C5_C8", "C6_C8", "C7_C8") ~ "less")) %>%
  group_by(class) %>% mutate(DE_Num=sum(DE_num), No_DE_Num=sum(No_DE_num)) %>% ungroup() %>%
  select(class, DE_Num, No_DE_Num) %>% as.data.frame() %>% distinct() %>% rowwise() %>%
  mutate(DE_freq = DE_Num / (DE_Num + No_DE_Num),
         No_DE_freq = No_DE_Num / (DE_Num + No_DE_Num)) %>% ungroup() %>% as.data.frame()
res1

chisq.test(matrix(c(622, 2656, 2194, 4362), nrow = 2)) # p-value < 2.2e-16

#= bootstrap
times=1000
res_bootstrap <- lapply(seq(times), function(t){
  #t <- 1
  sim_res <- res %>% sample_n(size = n(), replace = TRUE) %>% group_by(contrast) %>%
    mutate(DE_num=length(gene)) %>% ungroup() %>% select(contrast,DE_num) %>% distinct() %>% as.data.frame() %>%
    mutate(class = case_when(contrast %in% c("C5_C7", "C6_C7") ~ "more",
                             contrast %in% c("C5_C6", "C5_C8", "C6_C8", "C7_C8") ~ "less")) %>%
    group_by(class) %>% mutate(DE_Num=sum(DE_num)) %>% ungroup() %>% select(class, DE_Num) %>% 
    as.data.frame() %>% distinct() %>% mutate(DE_freq = ifelse(class=="more", DE_Num/3278, DE_Num/6556),
                                              times = t)
  return(sim_res)
}) %>% bind_rows()

res_bootstrap_new <- res_bootstrap %>% group_by(class) %>% 
  summarise(CI_lower1=quantile(DE_Num, 0.025),   # 2.5% 分位数
            CI_upper1=quantile(DE_Num, 0.975),   # 97.5% 分位数
            CI_lower2=quantile(DE_freq, 0.025),
            CI_upper2=quantile(DE_freq, 0.975),
            CI_method = "percentile") %>% as.data.frame()
res2 <- merge(res1, res_bootstrap_new, by="class") %>% mutate(CI_lower=DE_freq-CI_lower2,
                                                              CI_upper=CI_upper2-DE_freq)
res2$class <- factor(res2$class, levels = c("more", "less"))


pdf("/mnt/data5/disk/yangwj/Result_plots/Figure4_5/PASTRI/Figure5E.pdf",width=2.5,height=3.5)
ggplot(res2, aes(x=class, y=DE_freq))+geom_bar(stat = "identity",  position="identity", fill="#C6DBEF", width=0.6)+
  labs(x="Transition rate between a pair of cell states",
       y="Fraction of genes that are differentially expressed between the state pair")+ 
  geom_errorbar(aes(ymin=DE_freq-CI_lower, ymax=DE_freq+CI_upper),width=0.2, color="#4292C6")+
  theme_classic()+ scale_y_continuous(expand=c(0,0), limits = c(0, 0.36))+
  theme(axis.text.x=element_text(size=10,color="black",angle=30,hjust=1,vjust=1),
        axis.text.y=element_text(size=10,color="black"),axis.ticks=element_line(color="black"))
dev.off()





