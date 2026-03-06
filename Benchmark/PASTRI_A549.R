rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  library(PASTRI)
  library(parallel)
  library(ape)
  library(ggtree)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggmagnify)
  library(reshape2)
  library(ggbreak)
})

## ==================== source and define functions ===========================================
plot.TreeAndCelltype <- function(tree, tree_info, color_value, alphas){
  select.node.celltype <- tree_info %>% mutate(celltype=ifelse(cell_num>1, "multi_nodes", as.character(celltype))) %>%
    select(tiplabel, cell_num, celltype) %>% unique() %>% as.data.frame()
  rownames(select.node.celltype) <- select.node.celltype$tiplabel
  
  #==set use colors
  celltype.colors <- structure(c("#B2DF8A","#FE9999","#A4CEDE","#1F78B4","#3DA42D","#BBBFE0","gray"),
                               names=as.character(c(paste0("C",seq(0,5)), "multi_nodes")))
  #==only keep the celltype cols
  select.node.celltype.use <- as.data.frame(select.node.celltype[, "celltype", drop=F])
  select.node.celltype.use$celltype <- as.factor(select.node.celltype.use$celltype)
  
  # plot tree and celltypes
  p <- ggtree(tree,branch.length="none",layout="circular",size=0.1,color=color_value,alpha=alphas)+xlim(-10, NA)+ #circular or radial
    labs(title = "C6 (n = 2745)")
  #labs(title = "E8 (n = 1356)")
  p <- rotate_tree(p, 30)
  p <- 
    gheatmap(p, select.node.celltype.use, offset=0.05, width=0.1, colnames_angle=120, colnames_offset_y = .25, colnames = FALSE,color=NA) + xlim(-2, NA) +
    scale_fill_manual(name="cell types",values = celltype.colors)
  p
  return(p)
}


## ============ Run ===========================================================================
#== 0 plot tree
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/treeInfo.Rda")
tree_path <- "/mnt/data5/disk/wupeng/A549_lineage/2.PacBio_pepline/10.iqtree/filter_100indel/unique/"
C6.tree <- read.tree(paste0(tree_path, "C6/A549_C6", ".nwk"))
p.all.C6 <- plot.TreeAndCelltype(tree=C6.tree, tree_info=C6_treeInfo, color_value="black", alphas=0.3)

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/Figure7E.pdf",width=3,height=3)
p.all.C6
dev.off()


#== 1 get optimal transition
#calculate_lca_depths
tree_path <- "/mnt/data5/disk/wupeng/A549_lineage/2.PacBio_pepline/10.iqtree/filter_100indel/unique/"
C6.node.depth <- calculate_lca_depths(paste0(tree_path, "C6/A549_C6", ".nwk"))
E8.node.depth <- calculate_lca_depths(paste0(tree_path, "all_E8/A549_E8", ".nwk"))
G5.node.depth <- calculate_lca_depths(paste0(tree_path, "G5/A549_G5", ".nwk"))

save(C6.node.depth, E8.node.depth, G5.node.depth,
     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/TumorA549.Node.depth.Rda") 

#load some files
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/TumorA549.Node.depth.Rda")
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/treeInfo.Rda")

#get optimal trans
Sample <- c("C6", "E8", "G5") 
depth <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

df_mrca.mheight <- lapply(seq(length(Sample)), function(S){
  #S <- 1
  sel.Sample <- Sample[S]
  print(sel.Sample)
  sel.Onenode.depth <- paste0(sel.Sample,".node.depth") %>% get()
  ##get cell info
  sel.cell <- paste0(sel.Sample,"_treeInfo")  %>% get() %>% select(tiplabel, BC, celltype, cell_num)
  colnames(sel.cell) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
  
  ##
  cell_type <- c("C0","C1","C2","C3","C4","C5")
  Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
  
  df_depth <- lapply(seq(length(depth)), function(d){
    #d <- 4
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
      mutate(sample=sel.Sample)
    return(mrca.mheight_optimal_results_d1)
  }) %>% bind_rows()
  ##
  out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/"
  saveRDS(df_depth, file= paste0(out_path, sel.Sample,".op.trans",".Rds"))
  print("Successful")
  return(df_depth)
})


#== 2 compare correlation between samples
#rm(list=ls())
depth <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
C6.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/C6.op.trans.Rds")
E8.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/E8.op.trans.Rds")

cor.res <- lapply(seq(length(depth)), function(d){
  #d <- 2
  sel.depth <- depth[d]
  sel.depth.1 <- paste0("d_", sel.depth)
  print(sel.depth.1)
  ##
  sel.C6.op <- C6.op %>% filter(depth==sel.depth.1) %>% select(type_combined_new, start_cell, end_cell, norm_optimal_T)
  colnames(sel.C6.op) <- c("type.combined.new","start.cell","end.cell","C6.norm.optimal.T")
  sel.E8.op <- E8.op %>% filter(depth==sel.depth.1) %>% select(type_combined_new, norm_optimal_T)
  colnames(sel.E8.op) <- c("type.combined.new","E8.norm.optimal.T")
  
  all <- left_join(sel.C6.op, sel.E8.op, by="type.combined.new")
  ##compare
  C6_E8.cor <- cor.test(all$C6.norm.optimal.T, all$E8.norm.optimal.T, method = "pearson")
  C6_E8.rho <- cor.test(all$C6.norm.optimal.T, all$E8.norm.optimal.T, method = "spearman")
  C6_E8.fdist <- dist(rbind(all$C6.norm.optimal.T, all$E8.norm.optimal.T), method = "euclidean") %>% as.vector()
  C6_E8.mhdist <- dist(rbind(all$C6.norm.optimal.T, all$E8.norm.optimal.T), method = "manhattan") %>% as.vector()
  
  ##
  return(data.frame(depth=sel.depth,
                    times="real", 
                    stringsAsFactors =F,
                    C6_E8.cor=C6_E8.cor$estimate, C6_E8.cor.pval=C6_E8.cor$p.value,
                    C6_E8.rho=C6_E8.rho$estimate, C6_E8.rho.pval=C6_E8.rho$p.value,
                    C6_E8.fdist=C6_E8.fdist, C6_E8.mhdist=C6_E8.mhdist))
}) %>% bind_rows()

saveRDS(cor.res, "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/compare_real.Rds")

p.cor.res <- lapply(seq(length(depth)), function(d){
  #d <- 5
  sel.depth <- depth[d]
  sel.depth.1 <- paste0("d_", sel.depth)
  print(sel.depth.1)
  ##
  sel.C6.op <- C6.op %>% filter(depth==sel.depth.1) %>% select(type_combined_new, start_cell, end_cell, norm_optimal_T)
  colnames(sel.C6.op) <- c("type.combined.new","start.cell","end.cell","C6.norm.optimal.T")
  sel.E8.op <- E8.op %>% filter(depth==sel.depth.1) %>% select(type_combined_new, norm_optimal_T)
  colnames(sel.E8.op) <- c("type.combined.new","E8.norm.optimal.T")
  
  all <- left_join(sel.C6.op, sel.E8.op, by="type.combined.new")
  ##
  sel.cor.res <- cor.res %>% filter(depth==sel.depth)
  
  ##
  p.C6_E8 <- ggplot(all, aes(x = C6.norm.optimal.T, y = E8.norm.optimal.T))+
    geom_point(aes(color=start.cell, shape=end.cell), size=4)+
    labs(title=paste0(sel.depth.1, "[mheight]",
                      "\ncor = ", round(sel.cor.res$C6_E8.cor, 2), " P = ", round(log10(sel.cor.res$C6_E8.cor.pval),2),
                      "\neuclidean distance = ", round(sel.cor.res$C6_E8.fdist, 2)),
         x="PASTRI inferred transition rate of C6",y="PASTRI inferred transition rate of E8")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+
    xlim(0,0.6)+ylim(0,0.6)+scale_shape_manual(values=c(0,6,7,15,16,17))+
    scale_color_manual(values=c("#B2DF8A","#FE9999","#A4CEDE","#1F78B4","#3DA42D","#BBBFE0"))+
    theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
          legend.text=element_text(size=10,color="black"))
  return(p.C6_E8)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/PASTRI/Figure7G.pdf",width=4,height=5) 
p.cor.res[[8]]
dev.off()

#== 3 PASTRI accuracy
## 3.1 shuffle cellInfo and re-calculate transition rate by PASTRI
rm(list=ls())
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/TumorA549.Node.depth.Rda")
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig_TumorA549/treeInfo.Rda")
depth <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
sel.sample <- "C6"
other.sample.1 <- "E8"

n <- 1000
df_random_new1 <- lapply(seq(n), function(i){
  #i <- 1
  repeat{
    cat("times", i, "\n")
    sel.Onenode.depth <- paste0(sel.sample,".node.depth") %>% get() 
    other1.Onenode.depth <- paste0(other.sample.1,".node.depth") %>% get()
    
    ##load sel sample celltype infos and random tree leaf
    sel.celltype <- paste0(sel.sample, "_treeInfo") %>% get() %>% select(tiplabel, BC, celltype, cell_num)
    colnames(sel.celltype) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
    other1.celltype <- paste0(other.sample.1, "_treeInfo") %>% get() %>% select(tiplabel, BC, celltype, cell_num)
    colnames(other1.celltype) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
    
    ##random tree leaf
    sel.sim.One.celltype <- sel.celltype %>% filter(cellNum==1) %>% 
      mutate(celltype=sample(sel.celltype$celltype, size=nrow(.), replace=F))
    #table(sel.sim.One.celltype$celltype)
    other1.sim.One.celltype <- other1.celltype %>% filter(cellNum==1) %>% 
      mutate(celltype=sample(other1.celltype$celltype, size=nrow(.), replace=F))
    ##
    ##
    cell_type <- c("C0","C1","C2","C3","C4","C5")
    Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
    
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
        ##real
        C6.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/Re.C6.op.trans.Rds")
        E8.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/Re.E8.op.trans.Rds")
        
        sel_real_C6 <- C6.op %>% filter(depth==paste0("d_",sel.depth)) %>% mutate(sample="real_C6")
        sel_real_E8 <- E8.op %>% filter(depth==paste0("d_",sel.depth)) %>% mutate(sample="real_E8")
        ##compare correlation between samples
        ##
        #random
        sel.op <- sel_optimal_results_d %>% select(type_combined_new, start_cell, end_cell, norm_optimal_T)
        colnames(sel.op) <- c("type.combined.new","start.cell","end.cell", paste0(sel.sample, ".norm.optimal.T"))
        other_1.op <- other_1_optimal_results_d %>% select(type_combined_new, norm_optimal_T)
        colnames(other_1.op) <- c("type.combined.new", paste0(other.sample.1, ".norm.optimal.T"))
        #real
        C6.op <- sel_real_C6 %>% select(type_combined_new, norm_optimal_T)
        colnames(C6.op) <- c("type.combined.new", "real.C6.norm.optimal.T")
        E8.op <- sel_real_E8 %>% select(type_combined_new, norm_optimal_T)
        colnames(E8.op) <- c("type.combined.new", "real.E8.norm.optimal.T")
        #merge
        all <- full_join(full_join(sel.op, other_1.op, by="type.combined.new"),
                         full_join(C6.op, E8.op, by="type.combined.new"),by="type.combined.new")
        ##compare
        sim_C6_real_E8.cor <- cor.test(all$C6.norm.optimal.T, all$real.E8.norm.optimal.T, method = "pearson")
        sim_C6_real_E8.rho <- cor.test(all$C6.norm.optimal.T, all$real.E8.norm.optimal.T, method = "spearman")
        sim_C6_real_E8.fdist <- dist(rbind(all$C6.norm.optimal.T, all$real.E8.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_C6_real_E8.mhdist <- dist(rbind(all$C6.norm.optimal.T, all$real.E8.norm.optimal.T), method = "manhattan") %>% as.vector()
        
        sim_E8_real_C6.cor <- cor.test(all$E8.norm.optimal.T, all$real.C6.norm.optimal.T, method = "pearson")
        sim_E8_real_C6.rho <- cor.test(all$E8.norm.optimal.T, all$real.C6.norm.optimal.T, method = "spearman")
        sim_E8_real_C6.fdist <- dist(rbind(all$E8.norm.optimal.T, all$real.C6.norm.optimal.T), method = "euclidean") %>% as.vector()
        sim_E8_real_C6.mhdist <- dist(rbind(all$E8.norm.optimal.T, all$real.C6.norm.optimal.T), method = "manhattan") %>% as.vector()
        ##
        return(data.frame(depth=sel.depth,
                          times=i, 
                          stringsAsFactors =F,
                          sim_C6_real_E8.cor=sim_C6_real_E8.cor$estimate, sim_C6_real_E8.cor.pval=sim_C6_real_E8.cor$p.value,
                          sim_C6_real_E8.rho=sim_C6_real_E8.rho$estimate, sim_C6_real_E8.rho.pval=sim_C6_real_E8.rho$p.value,
                          sim_C6_real_E8.fdist=sim_C6_real_E8.fdist, sim_C6_real_E8.mhdist=sim_C6_real_E8.mhdist,
                          
                          sim_E8_real_C6.cor=sim_E8_real_C6.cor$estimate, sim_E8_real_C6.cor.pval=sim_E8_real_C6.cor$p.value,
                          sim_E8_real_C6.rho=sim_E8_real_C6.rho$estimate, sim_E8_real_C6.rho.pval=sim_E8_real_C6.rho$p.value,
                          sim_E8_real_C6.fdist=sim_E8_real_C6.fdist, sim_E8_real_C6.mhdist=sim_E8_real_C6.mhdist,
                          
                          C6_E8.cor=mean(c(sim_C6_real_E8.cor$estimate, sim_E8_real_C6.cor$estimate)),
                          C6_E8.rho=mean(c(sim_C6_real_E8.rho$estimate, sim_E8_real_C6.rho$estimate)),
                          C6_E8.fdist=mean(c(sim_C6_real_E8.fdist, sim_E8_real_C6.fdist)),
                          C6_E8.mhdist=mean(c(sim_C6_real_E8.mhdist, sim_E8_real_C6.mhdist))
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
out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/"
saveRDS(df_random_new1, file= paste0(out_path, "1000allDep_compare_shuffle_2mean_allLeafs_sim1.vs.real1",".Rds"))

## 3.2 compare with real transition rates (Figure 7)
rm(list=ls())
path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/"
real.res <- readRDS(paste0(path, "compare_real.Rds"))
sim.res <- readRDS(paste0(path, "1000allDep_compare_shuffle_2mean_allLeafs_sim1.vs.real1.Rds"))

depth <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

#get p value base on rank
sim.real.res_p_rank <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>% 
    select(depth, times, C6_E8.cor, C6_E8.rho, C6_E8.fdist, C6_E8.mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>%
    select(depth, times, C6_E8.cor, C6_E8.rho, C6_E8.fdist, C6_E8.mhdist)
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
  filter(methods=="C6_E8.fdist") %>% mutate(condition=paste0("<=",depth)) %>% select(-c(depth, methods))
colnames(sim.real.res_p_rank) <- c("p_rank", "condition")
sim.real.res_p_rank

#get p value base on test
sim.real.res_p.new <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(depth,C6_E8.cor,C6_E8.rho,C6_E8.fdist,C6_E8.mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(depth,C6_E8.cor,C6_E8.rho,C6_E8.fdist,C6_E8.mhdist)
  ##
  diff <- sel.op_real$C6_E8.fdist-mean(sel.op_shuffle$C6_E8.fdist)
  z_score <- (sel.op_real$C6_E8.fdist - mean(sel.op_shuffle$C6_E8.fdist)) / sd(sel.op_shuffle$C6_E8.fdist)
  
  ##
  cor.p <- wilcox.test(sel.op_shuffle$C6_E8.cor, mu=sel.op_real$C6_E8.cor)$p.value %>% -log10(.)
  fdist.p <- wilcox.test(sel.op_shuffle$C6_E8.fdist, mu=sel.op_real$C6_E8.fdist)$p.value %>% -log10(.)
  mhdist.p <- wilcox.test(sel.op_shuffle$C6_E8.mhdist, mu=sel.op_real$C6_E8.mhdist)$p.value %>% -log10(.)
  
  return(data.frame(depth = sel.depth, Diff = diff, z_score = z_score,
                    real.cor=sel.op_real$C6_E8.cor, cor.Pvalue=cor.p,
                    real.fdist=sel.op_real$C6_E8.fdist, fdist.Pvalue=fdist.p,
                    real.mhdist=sel.op_real$C6_E8.mhdist, mhdist.Pvalue=mhdist.p,
                    stringsAsFactors = F))
}) %>% bind_rows() %>%
  mutate(condition=paste0("<=",depth)) %>% #select(condition, Diff, z_score, real.fdist, fdist.Pvalue)
  select(condition, z_score, fdist.Pvalue, real.fdist)
##
sim.real.res_p.new <- full_join(sim.real.res_p.new, sim.real.res_p_rank, by="condition")
sim.real.res_p.new

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/PASTRI/Figure7F.new1.pdf",width=5.5,height=3.5)
ggplot(sim.real.res_p.new, aes(x=condition,y=real.fdist,color=condition))+geom_point(size=3)+#color="#1D91C0"
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Repeatability of robust\nAverage distance inferred from biological repeats")+
  theme_classic()+#scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.8", "#91C893", "#D6D6D6")) +
  scale_y_break(c(0.215,0.73),scales = 3/7, space=0.2)+
  #scale_y_continuous(limits = c(NA, max(sim.real.res_p.new$real.fdist, na.rm = TRUE) * 1.01))+
  theme(axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"),
        legend.position="none")
dev.off()


#====== Figure S8F PhyloVelo ==========
values <- 1.5535267824789065
pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/PASTRI/FigureS8F.pdf",width=6,height=5)
ggplot(sim.real.res_p.new, aes(x=condition,y=real.fdist,color=condition))+geom_point(size=6)+#color="#1D91C0"
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Repeatability of robust\ndistance inferred from biological repeats")+
  theme_classic()+geom_hline(yintercept=values, color="#E41A1C", linetype="dashed", linewidth=0.6) +
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.8", "#91C893", "#D6D6D6")) +
  scale_y_break(c(0.215,0.73),scales = 3/7, space=0.2)+
  #scale_y_continuous(limits = c(NA, max(sim.real.res_p.new$real.fdist, na.rm = TRUE) * 1.01))+
  theme(axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"),
        legend.position="none")
dev.off()


#====== Figure S8E significance test ==========
res.plot <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(depth, C6_E8.fdist)
  #select(depth, C6_E8.cor,C6_E8.rho,C6_E8.fdist,C6_E8.mhdist)
  op_shuffle.1 <- sel.op_shuffle %>% melt(., id="depth")
  
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(depth, C6_E8.fdist)
  #select(depth, C6_E8.cor,C6_E8.rho,C6_E8.fdist,C6_E8.mhdist)
  
  #set fact label name
  facet_labels <- sim.real.res_p.new %>% filter(depth==sel.depth) %>% 
    mutate(label=paste0(condition, " fdist ", " P = 10e-", round(fdist.Pvalue,2))) 
  facet_labels <- setNames(facet_labels$label, "C6_E8.fdist") 
  #plot
  p <- ggplot(op_shuffle.1, aes(x=value))+geom_histogram(binwidth=0.008,fill="#92C5DE")+
    geom_vline(data = sel.op_real, aes(xintercept = C6_E8.fdist), linetype = "dashed", color="#D73027")+
    facet_wrap(~ variable, scales = "free", labeller = as_labeller(facet_labels), nrow=1)+
    theme_classic()+scale_y_continuous(expand=c(0,0))+ labs(x="",y="Count")+ #theme_bw()+
    theme(axis.title=element_text(size=10,color="black"), axis.text=element_text(size=10,color="black"),
          axis.line=element_line(color="black"), strip.text=element_text(size=8,color="black"),
          #strip.background=element_rect(fill="#E6E6E6"),strip.text.x=element_text(margin=unit(rep(8,4),"pt")),
          panel.spacing = unit(1, "lines"))
  return(p)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/PASTRI/FigureS8E.pdf",width=3,height=3)
res.plot[[8]]
dev.off()


#====== Fig 7I ======
rm(list=ls())
C6.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/C6.op.trans.Rds") %>%
  filter(depth=="d_0.4") %>% select(-c(depth,sample))
colnames(C6.op) <- c("type_combined_new","start_cell","end_cell","C6")
E8.op <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig7_TumorA549/E8.op.trans.Rds") %>%
  filter(depth=="d_0.4") %>% select(type_combined_new, norm_optimal_T)
colnames(E8.op) <- c("type_combined_new","E8")

all <- full_join(C6.op, E8.op, by="type_combined_new") %>%
  mutate(mean_value=apply(.[,c("C6", "E8")], 1, mean) %>% round(., 4), 
         sd_value=apply(.[,c("C6", "E8")], 1, sd) %>% round(., 4)) %>% arrange(start_cell)

all %>% select(start_cell, end_cell, mean_value) %>% dcast(.,end_cell~start_cell,value.var="mean_value")

library(circlize)
net <- all %>% select(start_cell, end_cell, mean_value)
color <- structure(c("#B2DF8A","#FE9999","#A4CEDE","#1F78B4","#3DA42D","#BBBFE0"),names=c("C0","C1","C2","C3","C4","C5"))
group <- structure(as.character(rep(1:6)),names=c("C0","C1","C2","C3","C4","C5"))

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure7/PASTRI/Figure7I.pdf",width=3,height=3)
circos.clear()
chordDiagram(net, grid.col = color, directional=1, direction.type="diffHeight+arrows", 
             annotationTrack=c("name", "grid"), annotationTrackHeight=mm_h(c(1, 2)), link.sort=T,
             transparency=0.3, link.arr.lwd=net$mean_value*6, link.arr.length=0.2,
             #link.arr.col = "#B3B3B3", 
             link.arr.col = "#807F81", 
             scale=T, group = group, big.gap = 5)
dev.off()



