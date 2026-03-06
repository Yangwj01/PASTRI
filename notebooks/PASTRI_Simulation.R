rm(list = ls())
## ==================== import some librarys ==================================================
suppressMessages({
  #library(tidyverse)
  library(PASTRI)
  library(parallel)
  library(dplyr)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(tidyr)
  library(nortest)
  library(ggbreak)
  library(patchwork)
})

## ==================== source and define functions ===========================================
get_tans_mean <- function(trans, time){
  #trans <- trans_depth
  #time <- times
  trans_df <- mclapply(seq(time),function(t){
    #t <- 1
    trans_op_1 <- trans_depth[[t]]
    df_depth <- mclapply(seq(length(trans_op_1)), function(i){
      #i <- 1
      trans_op_2 <- trans_op_1[[i]]$optimal_norm_df_dataframe %>% mutate(sim_real=t)
      return(trans_op_2)
    }) %>% bind_rows()
    return(df_depth)
  }, mc.cores = 60) %>% bind_rows()
  
  ##mean
  dep <- unique(trans_df$depth)
  trans_df_mean <- mclapply(seq(length(dep)),function(d){
    #d <- 1
    sel.dep <- dep[d]
    sel.trans <- trans_df %>% filter(depth==sel.dep) %>% group_by(type_combined_new) %>%
      mutate(mean_norm_optimal_T=mean(norm_optimal_T), 
             sd_norm_optimal_T=sd(norm_optimal_T),
             se_norm_optimal_T = sd_norm_optimal_T / sqrt(n())) %>% 
      select(-c(norm_optimal_T, sim_real)) %>% distinct() %>% ungroup()
    return(sel.trans)
  }) %>% bind_rows()
  return(trans_df_mean)
}

## ============ Run ===========================================================================
###simulation 1000 trees and shuffle one times
#load some files
#type <- "ChainLike_ir"
type <- "ChainLike_r"
#type <- "Branched_ir"
#type <- "Branched_r"
#type <- "TwoPhase_ChainLike_ir"
#type <- "TwoPhase_ChainLike_r"
#type <- "TwoPhase_Branched_ir"
#type <- "TwoPhase_Branched_r"
times <- 1000

#== 1. Get the normalized depth of tip (cell) pairs
cat(type, "start", "\n")
node_depth <- mclapply(seq(times),function(t){
  #t <- 1
  sel.t <- t-1
  
  cat(Sys.time(), "time", sel.t, "\n")
  
  tree_path <- paste0("/mnt/data5/disk/yangwj/Scripts/Fig2/SmallTrees/PASTRI/",type,"/",type,"_",sel.t,"_tree.nwk")
  #Tree <- read.tree(tree_path)
  #ggtree(Tree)
  sample.node.depth <- calculate_lca_depths(tree_path)
  #head(sample.node.depth)
  #unique(sample.node.depth$lca_normalized_height)
  # output
  save(sample.node.depth, file=paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig2_simulation/SmallTrees/PASTRI_",
                                      type, "/", type,"_",sel.t,"_nodeDepth.Rda"))
  return(sample.node.depth)
},mc.cores = 5) 
cat(type, "success", "\n")


#== 2. Get the optimal transition matrix
trans_depth <- mclapply(seq(times),function(t){
  #t <- 1
  sel.t <- t-1
  cat(Sys.time(), "time", sel.t, "\n")
  #== 2.1 Get the optimal transition matrix
  cell_path <- paste0("/mnt/data5/disk/yangwj/Scripts/Fig2/SmallTrees/PASTRI/",type,"/",type,"_",sel.t,"_nodeInfos.txt")
  cell <- read.table(cell_path, header=T) %>% filter(grepl("^cell", nodeLabel))
  
  if (type %in% c("TwoPhase_ChainLike_ir", "TwoPhase_ChainLike_r", "TwoPhase_Branched_ir", "TwoPhase_Branched_r")) {
    depth <- c(1, 2)
  } else {
    depth <- c(1)
  }
  cell_type <- c("A", "B", "C", "D", "E")
  
  #==2.1 Set Bound_matrix
  if (type == "ChainLike_ir") {
    Bound_matrix <- matrix(c(Inf, Inf, 0, 0, 0,
                             0, Inf, Inf, 0, 0,
                             0, 0, Inf, Inf, 0,
                             0, 0, 0, Inf, Inf,
                             0, 0, 0, 0, Inf), 
                           nrow = length(cell_type), ncol = length(cell_type))
  } else if (type == "ChainLike_r") {
    Bound_matrix <- matrix(c(Inf, Inf, 0, 0, 0,
                             Inf, Inf, Inf, 0, 0,
                             0, Inf, Inf, Inf, 0,
                             0, 0, Inf, Inf, Inf,
                             0, 0, 0, Inf, Inf),
                           nrow = length(cell_type), ncol = length(cell_type))
  } else if (type == "Branched_ir") {
    Bound_matrix <- matrix(c(Inf, Inf, Inf, 0, 0,
                             0, Inf, 0, Inf, 0,
                             0, 0, Inf, 0, Inf,
                             0, 0, 0, Inf, 0,
                             0, 0, 0, 0, Inf),
                           nrow = length(cell_type), ncol = length(cell_type))
  } else if (type == "Branched_r") {
    Bound_matrix <- matrix(c(Inf, Inf, Inf, 0, 0,
                             Inf, Inf, 0, Inf, 0,
                             Inf, 0, Inf, 0, Inf,
                             0, Inf, 0, Inf, 0,
                             0, 0, Inf, 0, Inf),
                           nrow = length(cell_type), ncol = length(cell_type))
  } else {
    Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
  }
  #print(Bound_matrix)
  
  #==2.2 get optimal transition matrix per depth threshold
  nodeDepth_path <- paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig2_simulation/SmallTrees/PASTRI_",
                           type, "/", type,"_",sel.t,"_nodeDepth.Rda")
  load(nodeDepth_path)
  
  df_depth <- lapply(seq(length(depth)), function(d) {
    #d <- 1
    sel.depth <- depth[d]
    print(sel.depth)  # Print the current depth threshold
    
    # Call get_optimal_transition_matrix() with the specified parameters
    mrca.mheight_optimal_results_d <- get_optimal_transition_matrix(
      node_pair_depth = sample.node.depth,
      cell_info = cell,
      #Sel_u = "lca_normalized_height",
      Sel_u = "lca_depth",
      fi_depth = sel.depth,
      Bound_Matrix = Bound_matrix,
      cell_type_list = cell_type,
      mc.cores = 80
    )
    return(mrca.mheight_optimal_results_d)
  })
  return(df_depth)
})
trans_depth[[1]][[1]]$optimal_matrix
saveRDS(trans_depth, file=paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig2_simulation/SmallTrees/", type, ".lcaDep.optimal.trans.Rds"))


#== 3. Compare with actual transition matrix
#== 3.1 get actual transition rate
if (type == "ChainLike_ir") {
  T_true = matrix(c(0.73, 0.27, 0, 0, 0, 0, 0.85, 0.15, 0, 0, 0, 0, 0.90, 0.10, 0, 0, 0, 0, 0.62, 0.38, 0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else if (type == "ChainLike_r") {
  T_true = matrix(c(0.6, 0.4, 0, 0, 0, 0.2, 0.56, 0.24, 0, 0, 0, 0.12, 0.5, 0.38, 0, 0, 0, 0.3, 0.54, 0.16, 0, 0, 0, 0.37, 0.63), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else if (type == "Branched_ir") {
  T_true = matrix(c(0.5, 0.22, 0.28, 0, 0, 0, 0.71, 0, 0.29, 0, 0, 0, 0.86, 0, 0.14, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else if (type == "Branched_r") {
  T_true = matrix(c(0.5, 0.22, 0.28, 0, 0, 0.12, 0.59, 0, 0.29, 0, 0.18, 0, 0.68, 0, 0.14, 0, 0.21, 0, 0.79, 0, 0, 0, 0.17, 0, 0.83), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else if (type == "TwoPhase_ChainLike_ir") {
  # stage 2
  T_true = matrix(c(0.84, 0.16, 0, 0, 0, 0, 0.81, 0.19, 0, 0, 0, 0, 0.90, 0.10, 0, 0, 0, 0, 0.75, 0.25, 0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
  # stage 1
  #T_true = matrix(c(0.21, 0.16, 0.18, 0.13, 0.32, 0.28, 0.16, 0.21, 0.22, 0.13, 0.34, 0.06, 0.08, 0.22, 0.30, 0.20, 0.17, 0.14, 0.23, 0.26, 0.21, 0.15, 0.09, 0.25, 0.30), 
  #                nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else if (type == "TwoPhase_ChainLike_r") {
  # stage 2
  T_true = matrix(c(0.84, 0.16, 0, 0, 0, 0.08, 0.77, 0.15, 0, 0, 0, 0.20, 0.65, 0.15, 0, 0, 0, 0.12, 0.67, 0.21, 0, 0, 0, 0.16, 0.84), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
  # stage 1
  #T_true = matrix(c(0.21, 0.16, 0.18, 0.13, 0.32, 0.28, 0.16, 0.21, 0.22, 0.13, 0.34, 0.06, 0.08, 0.22, 0.30, 0.20, 0.17, 0.14, 0.23, 0.26, 0.21, 0.15, 0.09, 0.25, 0.30), 
  #                nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else if (type == "TwoPhase_Branched_ir") {
  # stage 2
  T_true = matrix(c(0.74, 0.16, 0.10, 0, 0, 0, 0.81, 0, 0.19, 0, 0, 0, 0.90, 0, 0.10, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
  # stage 1
  #T_true = matrix(c(0.21, 0.16, 0.18, 0.13, 0.32, 0.28, 0.16, 0.21, 0.22, 0.13, 0.34, 0.06, 0.08, 0.22, 0.30, 0.20, 0.17, 0.14, 0.23, 0.26, 0.21, 0.15, 0.09, 0.25, 0.30), 
  #                nrow = 5, ncol = 5, byrow = TRUE) %>% t()
} else {
  # stage 2
  T_true = matrix(c(0.74, 0.16, 0.10, 0, 0, 0.11, 0.70, 0, 0.19, 0, 0.13, 0, 0.77, 0, 0.10, 0, 0.15, 0, 0.85, 0, 0, 0, 0.08, 0, 0.92), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
  # stage 1
  T_true = matrix(c(0.21, 0.16, 0.18, 0.13, 0.32, 0.28, 0.16, 0.21, 0.22, 0.13, 0.34, 0.06, 0.08, 0.22, 0.30, 0.20, 0.17, 0.14, 0.23, 0.26, 0.21, 0.15, 0.09, 0.25, 0.30), 
                  nrow = 5, ncol = 5, byrow = TRUE) %>% t()
}

rownames(T_true) <- c("A","B","C","D","E")
colnames(T_true) <- c("A","B","C","D","E")
T_true
combinations <- expand.grid(start_cell=colnames(T_true), end_cell=rownames(T_true))
T_true <-mclapply(seq(nrow(combinations)),function(c){
  #c <- 1
  row_name <- combinations$end_cell[c]
  col_name <- combinations$start_cell[c]
  value <-T_true[row_name, col_name]  
  res <- data.frame(start_cell = col_name, end_cell = row_name, T_true=value) %>%
    mutate(type_combined_new=paste(start_cell, end_cell, sep = "_")) %>%
    select(type_combined_new,start_cell,end_cell,T_true)
  return(res)
}) %>% bind_rows()

#== 3.2 Compare
#=3.2.1 get pearson and Euclidean distance
trans_depth <- readRDS(paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig2_simulation/SmallTrees/", type, ".lcaDep.optimal.trans.Rds"))

cor.res <- mclapply(seq(times),function(t){
  #t <- 1
  cat(Sys.time(), "time", t, "\n")
  
  df_depth <- trans_depth[[t]]
  cor_res_df <- mclapply(seq(length(df_depth)), function(i){
    #i <- 1
    T_infer <- df_depth[[i]]$optimal_norm_df_dataframe
    True_infer <- left_join(T_true, T_infer, by=c("type_combined_new","start_cell", "end_cell")) %>%
      filter(T_true > 0)
    ##compare
    True_infer_cor <- cor.test(True_infer$T_true, True_infer$norm_optimal_T, method = "pearson")
    True_infer_fdist <- dist(rbind(True_infer$T_true, True_infer$norm_optimal_T), method = "euclidean") %>% as.vector()
    ##
    return(data.frame(depth=unique(T_infer$depth),
                      sim_tree=paste0("treeID", t-1),
                      times="real", 
                      stringsAsFactors =F,
                      True_infer.cor=True_infer_cor$estimate, 
                      True_infer.cor.pval=True_infer_cor$p.value,
                      True_infer.fdist=True_infer_fdist))
  }) %>% bind_rows()
  return(cor_res_df)
}, mc.cores = 60) %>% bind_rows()
cor.res %>% filter(depth=="d_1") %>% select(True_infer.fdist) %>% pull() %>% summary()

#=3.2.2 get mean Transition rate of all 1000 simulation trees per depth
tans_mean <- get_tans_mean(trans=trans_depth, time=times)

Dep <- unique(tans_mean$depth)

p.cor.res <- lapply(seq(length(Dep)), function(j){
  #j <- 1
  sel.dep <- Dep[j]
  T_infer <- tans_mean %>% filter(depth==sel.dep)
  True_infer <- left_join(T_true, T_infer, by=c("type_combined_new","start_cell", "end_cell")) %>%
    filter(T_true > 0)
  ##compare
  True_infer_cor.res <- cor.res %>% filter(depth==sel.dep)
  
  p.Ture_infer <- ggplot(True_infer, aes(x = T_true, y = mean_norm_optimal_T))+
    geom_point(aes(color=start_cell, shape=end_cell), size=3.8)+
    labs(title=paste0(True_infer$depth, 
                      "\ncor = ", round(mean(True_infer_cor.res$True_infer.cor), 2), " P=", log10(mean(True_infer_cor.res$True_infer.cor.pval)), 
                      "\nEuclidean distance = ", round(mean(True_infer_cor.res$True_infer.fdist), 2)),
         x="Actual transition rate", y="PASTRI inferred transition rate")+
    theme_classic()+geom_abline(intercept=0,slope=1,color="red",linetype="dashed")+coord_fixed()+
    geom_errorbar(aes(ymin=mean_norm_optimal_T-sd_norm_optimal_T, ymax=mean_norm_optimal_T+sd_norm_optimal_T),color="gray",width=0)+
    scale_shape_manual(values=c(0,6,15,16,17,18))+xlim(0,1)+ylim(0,1)+
    scale_color_manual(values = c("#84C3B8", "#7BA3C3", "#B8B3D1", "#EBA961", "#E1786C"))+
    theme(plot.title = element_text(size=10,color="black"),
          axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"))
  return(p.Ture_infer)
})

if (type %in% c("ChainLike_ir", "ChainLike_r", "Branched_ir", "Branched_r")) {
  pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/Figure2B_",type,"_1000trees.pdf"),width=8,height=4)
  p.cor.res[[1]]
} else {
  #pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/Figure2B_",type,"_1000trees.pdf"),width=12,height=4)
  pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/FigureS2_",type,"_Phase1_1000trees.pdf"),width=12,height=4)
  #pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/FigureS2_",type,"_PhaseMean_1000trees.pdf"),width=12,height=4)
  p.cor.res[[1]] | p.cor.res[[2]]
}
dev.off()


#===================================================================================================
#== 4 shuffle leaf states and random select three trees
sel_t <- sample(0:999, 3)
sim_times <- 1000

sim_true_trans_3 <- lapply(seq(length(sel_t)), function(sel){
  #sel <- 1
  sel.t <- sel_t[sel]
  cat("tree_ID", sel.t, "\n")
  
  sim_true_trans <- mclapply(seq(sim_times), function(t){
    #t <- 1
    cat("time", t, "\n")
    #== 4.1.1 shuffle leaf states
    cell_path <- paste0("/mnt/data5/disk/yangwj/Scripts/Fig2/SmallTrees/PASTRI/",type,"/",type,"_",sel.t,"_nodeInfos.txt")
    cell <- read.table(cell_path, header=T) %>% filter(grepl("^cell", nodeLabel))
    cell_suffle <- cell %>% mutate(raw_celltype=celltype) %>% 
      mutate(celltype=sample(raw_celltype, size=nrow(.), replace=F))
    if (type %in% c("TwoPhase_ChainLike_ir", "TwoPhase_ChainLike_r", "TwoPhase_Branched_ir", "TwoPhase_Branched_r")) {
      depth <- c(1, 2)
    } else {
      depth <- c(1)
    }
    cell_type <- c("A", "B", "C", "D", "E")
    
    #==4.1.2 Set Bound_matrix
    if (type == "ChainLike_ir") {
      Bound_matrix <- matrix(c(Inf, Inf, 0, 0, 0,
                               0, Inf, Inf, 0, 0,
                               0, 0, Inf, Inf, 0,
                               0, 0, 0, Inf, Inf,
                               0, 0, 0, 0, Inf), 
                             nrow = length(cell_type), ncol = length(cell_type))
    } else if (type == "ChainLike_r") {
      Bound_matrix <- matrix(c(Inf, Inf, 0, 0, 0,
                               Inf, Inf, Inf, 0, 0,
                               0, Inf, Inf, Inf, 0,
                               0, 0, Inf, Inf, Inf,
                               0, 0, 0, Inf, Inf),
                             nrow = length(cell_type), ncol = length(cell_type))
    } else if (type == "Branched_ir") {
      Bound_matrix <- matrix(c(Inf, Inf, Inf, 0, 0,
                               0, Inf, 0, Inf, 0,
                               0, 0, Inf, 0, Inf,
                               0, 0, 0, Inf, 0,
                               0, 0, 0, 0, Inf),
                             nrow = length(cell_type), ncol = length(cell_type))
    } else if (type == "Branched_r") {
      Bound_matrix <- matrix(c(Inf, Inf, Inf, 0, 0,
                               Inf, Inf, 0, Inf, 0,
                               Inf, 0, Inf, 0, Inf,
                               0, Inf, 0, Inf, 0,
                               0, 0, Inf, 0, Inf),
                             nrow = length(cell_type), ncol = length(cell_type))
    } else {
      Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
    }
    #print(Bound_matrix)
    
    #==4.1.3 get optimal transition matrix per depth threshold
    nodeDepth_path <- paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig2_simulation/SmallTrees/PASTRI_",
                             type, "/", type,"_",sel.t,"_nodeDepth.Rda")
    load(nodeDepth_path)
    
    #==4.1.4 infer sim trsansition matrix
    sim_df_depth <- lapply(seq(length(depth)), function(d) {
      #d <- 1
      sel.depth <- depth[d]
      print(sel.depth)  # Print the current depth threshold
      
      # Call get_optimal_transition_matrix() with the specified parameters
      sim_mrca.mheight_optimal_results_d <- get_optimal_transition_matrix(
        node_pair_depth = sample.node.depth,
        cell_info = cell_suffle,
        #Sel_u = "lca_normalized_height",
        Sel_u = "lca_depth",
        fi_depth = sel.depth,
        Bound_Matrix = Bound_matrix,
        cell_type_list = cell_type,
        mc.cores = 80
      )
      return(sim_mrca.mheight_optimal_results_d)
    })
    
    #==4.1.5 compare sim trsansition matrix with T_true
    Ture_sim_cor.res <- lapply(seq(length(sim_df_depth)), function(s){
      #s <- 1
      sim_T_infer <- sim_df_depth[[s]]$optimal_norm_df_dataframe
      sim_True_infer <- left_join(T_true, sim_T_infer, by=c("type_combined_new","start_cell", "end_cell")) %>%
        filter(T_true > 0)
      ##compare
      True_infer_cor <- cor.test(sim_True_infer$T_true, sim_True_infer$norm_optimal_T, method = "pearson")
      True_infer_fdist <- dist(rbind(sim_True_infer$T_true, sim_True_infer$norm_optimal_T), method = "euclidean") %>% as.vector()
      ##
      return(data.frame(depth=depth[s],
                        sim_tree=paste0("treeID",sel.t ),
                        times=paste0("sim_", t), 
                        stringsAsFactors =F,
                        sim_infer.cor=True_infer_cor$estimate, 
                        sim_infer.cor.pval=True_infer_cor$p.value,
                        sim_infer.fdist=True_infer_fdist))
    }) %>% bind_rows()
    
    return(Ture_sim_cor.res)
  },mc.cores = 60) %>% bind_rows()
  
  return(sim_true_trans)
}) %>% bind_rows()

save(sim_true_trans_3, file=paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig2_simulation/SmallTrees/", type, ".shuffle.res.Rda"))

#== new 5 for random select three trees
#get p value base on rank
sim.real.res_p_rank <- lapply(seq(length(sel_t)), function(sel){
  #sel <- 1
  sel.t <- sel_t[sel]
  cat("tree_ID", sel.t, "\n")
  
  sim.real.res_p_rank1 <- lapply(seq(length(Dep)), function(d){
    #d <- 1
    sel.depth <- Dep[d]
    print(sel.depth)
    #
    sel.op_shuffle <- sim_true_trans_3 %>% mutate(depth=paste0("d_",depth)) %>% filter(depth==sel.depth) %>% 
      select(which(!grepl("pval", names(.)))) %>% filter(sim_tree==paste0("treeID",sel.t))
    colnames(sel.op_shuffle) <- c("depth","sim_tree" ,"times", "infer.cor", "infer.fdist")
    sel.op_real <- cor.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>%
      filter(sim_tree==paste0("treeID",sel.t)) 
    colnames(sel.op_real) <- c("depth","sim_tree" ,"times", "infer.cor", "infer.fdist")
    
    sel.real.sim <- rbind(sel.op_shuffle, sel.op_real)
    col_name <- colnames(sel.real.sim)  
    #
    res_p <- mclapply(4:length(col_name), function(col){
      #col <- 4
      sel.col <- col_name[col]
      if (grepl("cor|rho", sel.col)) {
        sel.data <- sel.real.sim %>% select(c(1:3,sym(sel.col))) %>% .[order(.[,4], decreasing = T),]
      } else {
        sel.data <- sel.real.sim %>% select(c(1:3,sym(sel.col))) %>% .[order(.[,4], decreasing = F),]
      }
      p.value <- which(grepl("real", sel.data$times))/1000 #; 100-0.01
      return(data.frame(depth = sel.depth,
                        sim_tree = paste0("treeID",sel.t),
                        methods=sel.col,
                        p_value=p.value,
                        stringsAsFactors = F))
    }) %>% bind_rows()
    return(res_p)
  }) %>% bind_rows() %>% 
    mutate(condition=paste0("<=",depth))
  return(sim.real.res_p_rank1)
}) %>% bind_rows() %>%
  filter(methods=="infer.fdist")
sel.cor.res <- cor.res %>% filter(sim_tree %in% unique(sim_true_trans_3$sim_tree)) %>% 
  select(depth, sim_tree, True_infer.fdist)
left_join(sim.real.res_p_rank, sel.cor.res, by=c("depth","sim_tree")) %>% filter(depth %in% c("d_1", "d_2"))

#plot
res.plot <- lapply(seq(length(Dep)), function(d){
  #d <- 2
  sel.depth <- Dep[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim_true_trans_3 %>% mutate(depth=paste0("d_",depth)) %>% filter(depth==sel.depth) %>%
    select(depth, sim_tree, sim_infer.fdist)
  sel.op_real <- cor.res %>% filter(depth==sel.depth & sim_tree %in% unique(sel.op_shuffle$sim_tree)) %>% 
    select(depth, sim_tree, True_infer.fdist)
  
  sel.P <- sim.real.res_p_rank %>% filter(depth==sel.depth)
  #plot
  max_density <- max(ggplot_build(ggplot(sel.op_shuffle, aes(x = sim_infer.fdist, fill = sim_tree)) + geom_density(alpha = 0.2, color = "gray20"))$data[[1]]$density)
  #cutoff_1 <- max(sel.op_real$True_infer.fdist) + 0.05 
  #cutoff_2 <- min(sel.op_shuffle$sim_infer.fdist)
  
  #"#2F83AF" "#3F8EAA" "#509AA6" "#60A5A1" "#70B19C" "#81BC98" "#91C893" "#A1D38E" "#B2DF8A"
  #"#F8766D","#00BA38","#619CFF"
  colors <- c("#00BA38", "#F8766D", "#1D91C0")
  res.p <- 
    ggplot(sel.op_shuffle, aes(x=sim_infer.fdist, fill=sim_tree))+ geom_density(alpha=0.2,color="gray20")+
    #geom_histogram(binwidth=0.005)+
    geom_segment(data=sel.op_real, 
                 aes(x=True_infer.fdist, xend=True_infer.fdist ,y=max_density*0.6, yend=0, color=sim_tree),alpha=0.8,
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"))+
    labs(title = paste0("Adjacencies at ", sel.depth, " P = ", sel.P$p_value),x="",y="Density")+
    #xlim(min(sel.op_real$True_infer.fdist)-0.05, max(sel.op_shuffle$sim_infer.fdist)+0.01)+
    #scale_x_break(c(cutoff_1,cutoff_2),scales = 3,space=0.2)+
    theme_classic()+scale_y_continuous(limits=c(0, max_density*1.1), expand=c(0,0))+
    #scale_fill_manual(values = c("#3F8EAA","pink","#A1D38E"))
    #scale_fill_manual(values = c("#1D91C0","#B2DF8A","#E5C494"))
    scale_fill_manual(values = colors)+ scale_color_manual(values =colors)+
    theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
          axis.text=element_text(size=10,color="black"), legend.position = "none",
          axis.line=element_line(color="black"), panel.spacing = unit(1, "lines"),
          axis.line.x.top=element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top=element_blank())
  return(res.p)
})

res.plot.new <- lapply(seq(length(Dep)), function(d){
  #d <- 2
  sel.depth <- Dep[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim_true_trans_3 %>% mutate(depth=paste0("d_",depth)) %>% filter(depth==sel.depth) %>%
    select(depth, sim_tree, sim_infer.fdist)
  sel.op_real <- cor.res %>% filter(depth==sel.depth & sim_tree %in% unique(sel.op_shuffle$sim_tree)) %>% 
    select(depth, sim_tree, True_infer.fdist)
  
  sel.P <- sim.real.res_p_rank %>% filter(depth==sel.depth)
  #plot
  max_density <- max(ggplot_build(ggplot(sel.op_shuffle, aes(x = sim_infer.fdist, fill = sim_tree)) + geom_density(alpha = 0.2, color = "gray20"))$data[[1]]$density)
  cutoff_1 <- max(sel.op_real$True_infer.fdist) + 0.05 
  cutoff_2 <- min(sel.op_shuffle$sim_infer.fdist)
  
  #"#2F83AF" "#3F8EAA" "#509AA6" "#60A5A1" "#70B19C" "#81BC98" "#91C893" "#A1D38E" "#B2DF8A"
  #"#F8766D","#00BA38","#619CFF"
  colors <- c("#00BA38", "#F8766D", "#1D91C0")
  res.p <- 
    ggplot(sel.op_shuffle, aes(x=sim_infer.fdist, fill=sim_tree))+ geom_density(alpha=0.2,color="gray20")+
    #geom_histogram(binwidth=0.005)+
    geom_segment(data=sel.op_real, 
                 aes(x=True_infer.fdist, xend=True_infer.fdist ,y=max_density*0.6, yend=0, color=sim_tree),alpha=0.8,
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"))+
    labs(title = paste0("Adjacencies at ", sel.depth, " P = ", sel.P$p_value),x="",y="Density")+
    xlim(min(sel.op_real$True_infer.fdist)-0.05, max(sel.op_shuffle$sim_infer.fdist)+0.01)+
    scale_x_break(c(cutoff_1,cutoff_2),scales = 3,space=0.2)+
    theme_classic()+scale_y_continuous(limits=c(0, max_density*1.1), expand=c(0,0))+
    #scale_fill_manual(values = c("#3F8EAA","pink","#A1D38E"))
    #scale_fill_manual(values = c("#1D91C0","#B2DF8A","#E5C494"))
    scale_fill_manual(values = colors)+ scale_color_manual(values =colors)+
    theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
          axis.text=element_text(size=10,color="black"), legend.position = "none",
          axis.line=element_line(color="black"), panel.spacing = unit(1, "lines"),
          axis.line.x.top=element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top=element_blank())
  return(res.p)
})

res.plot.phase1 <- lapply(seq(length(Dep)), function(d){
  #d <- 3
  sel.depth <- Dep[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim_true_trans_3 %>% mutate(depth=paste0("d_",depth)) %>% filter(depth==sel.depth) %>%
    select(depth, sim_tree, sim_infer.fdist)
  sel.op_real <- cor.res %>% filter(depth==sel.depth & sim_tree %in% unique(sel.op_shuffle$sim_tree)) %>% 
    select(depth, sim_tree, True_infer.fdist)
  
  sel.P <- sim.real.res_p_rank %>% filter(depth==sel.depth)
  #plot
  max_density <- max(ggplot_build(ggplot(sel.op_shuffle, aes(x = sim_infer.fdist, fill = sim_tree)) + geom_density(alpha = 0.2, color = "gray20"))$data[[1]]$density)
  cutoff_1 <- max(sel.op_shuffle$sim_infer.fdist)+0.01
  cutoff_2 <- min(sel.op_real$True_infer.fdist)- 0.05 
  
  colors <- c("#F8766D","#1D91C0","#00BA38")
  if(d==3){
    res.p <- 
      ggplot(sel.op_shuffle, aes(x=sim_infer.fdist, fill=sim_tree))+ geom_density(alpha=0.2,color="gray20")+
      #geom_histogram(binwidth=0.005)+
      geom_segment(data=sel.op_real, 
                   aes(x=True_infer.fdist, xend=True_infer.fdist ,y=max_density*0.6, yend=0, color=sim_tree),alpha=0.8,
                   arrow = arrow(length = unit(0.1, "cm"), type = "closed"))+
      labs(title = paste0("Adjacencies at D <= ", sel.depth, "P = ", sel.P$p_value),x="",y="Density")+
      xlim(min(sel.op_shuffle$sim_infer.fdist), max(sel.op_real$True_infer.fdist)+0.05)+
      theme_classic()+scale_y_continuous(limits=c(0, max_density*1.1), expand=c(0,0))+
      scale_fill_manual(values = colors)+ scale_color_manual(values =colors)+
      theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
            axis.text=element_text(size=10,color="black"), legend.position = "none",
            axis.line=element_line(color="black"), panel.spacing = unit(1, "lines"),
            axis.line.x.top=element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top=element_blank())
  } else {
    res.p <- 
      ggplot(sel.op_shuffle, aes(x=sim_infer.fdist, fill=sim_tree))+ geom_density(alpha=0.2,color="gray20")+
      #geom_histogram(binwidth=0.005)+
      geom_segment(data=sel.op_real, 
                   aes(x=True_infer.fdist, xend=True_infer.fdist ,y=max_density*0.6, yend=0, color=sim_tree),alpha=0.8,
                   arrow = arrow(length = unit(0.1, "cm"), type = "closed"))+
      labs(title = paste0("Adjacencies at D <= ", sel.depth, "P = ", sel.P$p_value),x="",y="Density")+
      xlim(min(sel.op_shuffle$sim_infer.fdist)-0.01, max(sel.op_real$True_infer.fdist)+0.05)+
      scale_x_break(c(cutoff_1,cutoff_2),scales = 0.2,space=0.2)+
      theme_classic()+scale_y_continuous(limits=c(0, max_density*1.1), expand=c(0,0))+
      #scale_fill_manual(values = c("#3F8EAA","pink","#A1D38E"))
      #scale_fill_manual(values = c("#1D91C0","#B2DF8A","#E5C494"))
      scale_fill_manual(values = colors)+ scale_color_manual(values =colors)+
      theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
            axis.text=element_text(size=10,color="black"), legend.position = "none",
            axis.line=element_line(color="black"), panel.spacing = unit(1, "lines"),
            axis.line.x.top=element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top=element_blank())
  }
  
  return(res.p)
})

if (type %in% c("ChainLike_ir", "ChainLike_r", "Branched_ir", "Branched_r")) {
  pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/Figure2D_",type,"_3trees.pdf"),width=2.8,height=3.08)
  res.plot[[1]]
} else {
  #pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/Figure2D_",type,"_3trees.pdf"),width=5.6,height=3.08)
  #res.plot.new[[1]]|res.plot.new[[2]]
  pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/FigureS2_",type,"_Phase1_3trees.pdf"),width=5.6,height=3.08)
  #pdf(paste0("/mnt/data5/disk/yangwj/Result_plots/Figure2/PASTRI/FigureS2_",type,"_PhaseMean_3trees.pdf"),width=5.6,height=6.16)
  res.plot.phase1[[1]]|res.plot.phase1[[2]]
  #res.plot.new[[1]]|res.plot.new[[2]]
}
dev.off()





