rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
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
treeDist.pergene <- function(file_infos){
  #file_infos=filter_file
  gene_list <- unique(file_infos$gene)
  res.tree.pergene <- mclapply(seq(length(gene_list)),function(g){
    #g <- 3
    sel.gene <- gene_list[g]
    sel.fileInfo <- file_infos %>% filter(gene==sel.gene)
    file_list_pergene <- unique(sel.fileInfo$filename)
    res.sel.gene <- mclapply(seq(length(file_list_pergene)),function(f){
      #f <-1
      sel.file <- file_list_pergene[f]
      #data <- read.csv(sel.file,stringsAsFactors=FALSE)
      #
      tree_path <- file.path("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/treefile_Celegans", paste0(sel.file, ".nwk"))
      #terminal.cell <- data  %>% select(cell,time) %>% filter(time==max(time))
      #== get tree dist
      onecell.node.depth <- calculate_lca_depths(tree_path)
      
      #== set precision to 8
      onecell.node.pair <- onecell.node.depth %>% mutate(lca_depth_norm.raw=lca_depth_norm) %>% 
        mutate(lca_depth_norm=round(lca_depth_norm,8)) #%>% group_by(lca_depth_norm) %>%
      #add_count(lca_depth_norm, name="num.cell.perDepth") %>% ungroup()
      onecell.node.pair <- onecell.node.pair %>% mutate(gene=sel.gene,file=sel.file)
      
      #== output treeDist file
      output_filename <- basename(sel.file)
      output_path <- file.path("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/treeDist", output_filename)
      write.csv(onecell.node.pair, file=output_path, row.names=F)
      return(data.frame(treeDist_path=output_path,gene=sel.gene,file=sel.file))
    }) %>% bind_rows()
    return(res.sel.gene)
  }) %>% bind_rows()
  return(res.tree.pergene)
}

get.optimal.irreversible.transT.allgene <- function(file_infos, all_treeDist, use.depth, sel.depth){
  #file_infos = filter_file
  #all_treeDist = all.treeDist
  #use.depth = "lca_depth_norm"
  #sel.depth = 0.2
  gene_list <- unique(file_infos$gene)#[26]
  optimal.trans.pergene <- mclapply(seq(length(gene_list)),function(g){
    #g <- 1
    #sel.gene <- "elt-7"
    sel.gene <- gene_list[g]
    sel.fileInfo <- file_infos %>% filter(gene==sel.gene)
    file_list_pergene <- unique(sel.fileInfo$filename)
    # Optimal trans probability
    optimal.trans.sel.gene <- mclapply(seq(length(file_list_pergene)),function(f){
      #f <- 1
      sel.file <- file_list_pergene[f]
      data <- read.csv(paste0("/mnt/data5/disk/yangwj/Scripts/Cell_State_Transition/Celegans_EPIC/",sel.file),stringsAsFactors=FALSE)
      #calculate trans probability
      #==1 get tree Dist file
      sel.treeDist.file.path <- all_treeDist %>% filter(file==sel.file) %>% select(treeDist_path) %>% pull()
      treeDist <- read.csv(sel.treeDist.file.path,stringsAsFactors=FALSE)
      ##
      treeDist <- treeDist[,1:6]
      colnames(treeDist) <- c("node1", "node2", "subtree_root", "lca_depth", "lca_depth_norm", "lca_normalized_height")
      
      #==2 annotation cell state base expression or Fluorescence intensity
      terminal.cell <- data %>% select("cell","time","blot") %>% group_by(cell) %>% 
        mutate(start.time=min(time),end.time=max(time)) %>% ungroup() %>% 
        #find the cells when time==max(time), and use the fluorescence at that time as the standard
        filter(time==max(time)) %>% select(-time) %>% 
        group_by(cell) %>% mutate(mean.blot.raw=mean(blot)) %>% ungroup() %>% 
        mutate(mean.blot=mean.blot.raw-min(mean.blot.raw)) %>% select(-blot) %>% unique() %>%
        mutate(state=case_when(mean.blot>(max(mean.blot)*0.3)~"ON",mean.blot<=(max(mean.blot)*0.3)~"OFF"))
      #table(terminal.cell$state)
      terminal.cell <- terminal.cell %>% mutate(cell_ID=cell, cellNum=1) %>% select(cell, cell_ID, state, cellNum)
      colnames(terminal.cell) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
      ##
      cell_type <- c("OFF","ON")
      Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
      Bound_matrix[upper.tri(Bound_matrix)] <- 0
      
      #==3 calculate_transMatrix
      optimal.transMatrix <- get_optimal_transition_matrix(
        node_pair_depth = treeDist,
        cell_info = terminal.cell,
        Sel_u = use.depth,
        fi_depth = sel.depth,
        Bound_Matrix = Bound_matrix,
        cell_type_list = cell_type,
        mc.cores = 80
      )
      
      optimal_trans <- optimal.transMatrix$optimal_norm_df_dataframe %>% 
        mutate(gene=sel.gene,
               file=sel.file,
               rep=paste("rep", f, sep=""),
               max.end.time=unique(terminal.cell$end.time)) 
      
      freq.state <- terminal.cell %>% mutate(celltype = factor(celltype, levels = c("OFF","ON"))) %>%
        `[`("celltype") %>% table() %>% as.data.frame() %>% mutate(freq = Freq/sum(Freq))
      colnames(freq.state) <- c("start_cell","start.num","start.freq")
      optimal.norm.T <- left_join(optimal_trans, freq.state, by="start_cell") %>% as.data.frame()
      return(optimal.norm.T)
    }) %>% bind_rows()
    return(optimal.trans.sel.gene)
  }) %>% bind_rows()
  return(optimal.trans.pergene)
}

get.sheer.trans.branch <- function(file_infos, use.depth, sel.depth){
  #file_infos = filter_file
  #use.depth = "lca_depth_norm"
  #sel.depth = 0.2
  gene_list <- unique(file_infos$gene)
  sheer.T.pergene <- mclapply(seq(length(gene_list)),function(g){
    #g <- 3
    #sel.gene="elt-7"
    #sel.gene="mep-1"
    #sel.gene <- "F38C2.7"
    sel.gene <- gene_list[g]
    sel.fileinfo <- file_infos %>% filter(gene==sel.gene)
    file_list <- unique(sel.fileinfo$filename)
    sheer.T.perfile <- mclapply(seq(length(file_list)),function(f){
      #f <- 1
      file_path <- file_list[f]
      #file_path <- "CD20070310_elt7_c2a.csv"
      #file_path <- "CD20080412_ama-1_3A3_5.csv"
      path <- "/mnt/data5/disk/yangwj/Scripts/Cell_State_Transition/Celegans_EPIC/"
      data <- read.csv(paste0(path, file_path)) %>% select("cell","time","blot") %>% group_by(cell) %>% 
        mutate(s.time=min(time),e.time=max(time)) %>% ungroup() %>% filter(time==e.time) %>% select(-time) %>%
        group_by(cell) %>% mutate(mean.blot.raw=mean(blot)) %>% ungroup() %>% 
        mutate(blot.new=mean.blot.raw-min(mean.blot.raw)) %>% select(-blot) %>% unique() %>%
        mutate(state=case_when(blot.new>(max(blot.new)*0.3)~"ON", blot.new<=(max(blot.new)*0.3)~"OFF"))
      #table(data$state)
      
      #==1 read tree nwk file and get two-cell (branch) trans
      tree_path <- file.path("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/treefile_Celegans", paste0(file_path, ".nwk"))
      tree <- read.tree(tree_path)
      df_tree <- fortify(tree)
      df_tree.new <- df_tree %>% mutate(parent.label=df_tree$label[match(.$parent, df_tree$node)]) %>%
        select(parent, node, parent.label, label) %>% filter(label!="P0")
      colnames(df_tree.new) <- c("parent","node","parent.label","node.label")
      df_tree.new <- df_tree.new %>% mutate(state.parent=data$state[match(.$parent.label, data$cell)],
                                            state.node=data$state[match(.$node.label, data$cell)],
                                            #e.time.parent=data$e.time[match(.$parent.label, data$cell)],
                                            #e.time.node=data$e.time[match(.$node.label, data$cell)]
      ) %>%
        mutate(type.combined=paste(.$state.parent, .$state.node, sep = "_")) %>%
        filter(state.parent!="NA" & state.node!="NA")
      #table(df_tree.new$type.combined)
      
      #==2 add distance u 
      all_parents <- unique(df_tree$parent)
      
      subtree.depth <- lapply(1:length(all_parents),function(j){
        current_parent <- all_parents[j]
        ## calculate Dr (distance to root)
        i <-  which(df_tree$node == current_parent)
        a <- 0
        repeat{
          if(df_tree$node[i]==df_tree$parent[i])
            break 
          else{x <- df_tree$node[i]; y <- df_tree$parent[i]; i <- which(df_tree$node == y); a=a+1;}} 
        ## calculate Ds (distance to  farthest leaf)
        daughter_df <- df_tree %>% dplyr::filter(parent %in% current_parent & (parent != node))
        b <- 0
        repeat{
          if(nrow(daughter_df) == 0) break 
          else{
            daughter_df <- df_tree %>% dplyr::filter((parent %in% daughter_df$node) & (parent != node)); b=b+1;}}
        ## return results 
        return(data.frame("subtree.root" = current_parent, "Dr" = a, "Ds" = b))
      }) %>% bind_rows()
      
      subtree.depth$Dt <- max(subtree.depth$Ds)
      subtree.depth$Ds.norm <- subtree.depth$Ds/subtree.depth$Dt
      subtree.depth$mdepth <- ((subtree.depth$Dr/subtree.depth$Dt)+(1-subtree.depth$Ds/subtree.depth$Dt))/2
      subtree.depth$node.height <- 1-subtree.depth$mdepth
      
      treeDe <- subtree.depth %>% select(subtree.root, Ds, Ds.norm, node.height)
      colnames(treeDe) <- c("subtree_root", "lca_depth", "lca_depth_norm", "lca_normalized_height")
      ##
      if (use.depth == "lca_normalized_height") {
        treeDe <- treeDe %>% mutate(u=lca_normalized_height)
      } else {
        treeDe <- treeDe %>% mutate(u=lca_depth_norm)
      }
      
      df_tree.new1 <- df_tree.new %>% mutate(D.u=treeDe$u[match(.$parent, treeDe$subtree_root)])
      
      #==3 calculate sheer trans
      sheer.trans <- df_tree.new1 %>% filter(D.u <= sel.depth) %>% group_by(type.combined) %>%
        add_count(type.combined, name="num") %>% ungroup() %>% select(state.parent,type.combined,num) %>% unique() %>%
        mutate(freq=num/sum(num)) %>% group_by(state.parent) %>% mutate(sheer.T=freq/sum(freq)) %>% ungroup()
      #mutate(sheer.T=freq) 
      sheer.trans <- sheer.trans %>% as.data.frame() %>% mutate(gene=sel.gene,file=file_path)
      return(sheer.trans)
    },mc.cores = 60) %>% bind_rows()
    return(sheer.T.perfile)
  }) %>% bind_rows()
  return(sheer.T.pergene)
}

get.sim.op.irreversible.transT.allgene <- function(file_infos, all_treeDist, use.depth, depth, n){
  #file_infos = filter_file
  #all_treeDist = all.treeDist
  #use.depth = "lca_depth_norm"
  #depth = c(0.2, 0.4, 0.6, 0.8, 1.0)
  #n = 1
  gene_list <- unique(file_infos$gene)
  sim.op.trans.pergene <- lapply(seq(length(gene_list)),function(g){
    #g <- 1
    #sel.gene <- "elt-7"
    #sel.gene <- "F38C2.7"
    sel.gene <- gene_list[g]
    sel.fileInfo <- file_infos %>% filter(gene==sel.gene)
    file_list_pergene <- unique(sel.fileInfo$filename)
    # Optimal trans probability
    sim.op.trans.sel.gene <- lapply(seq(length(file_list_pergene)),function(f){
      #f <- 1
      sel.file <- file_list_pergene[f]
      cat("gene", g, sel.gene, "file", f, sel.file, "\n")
      
      data <- read.csv(paste0("/mnt/data5/disk/yangwj/Scripts/Cell_State_Transition/Celegans_EPIC/", sel.file),stringsAsFactors=FALSE)
      #calculate trans probability
      #==1 get tree Dist file
      sel.treeDist.file.path <- all_treeDist %>% filter(file==sel.file) %>% select(treeDist_path) %>% pull()
      treeDist <- read.csv(sel.treeDist.file.path,stringsAsFactors=FALSE)
      ##
      treeDist <- treeDist[,1:6]
      colnames(treeDist) <- c("node1", "node2", "subtree_root", "lca_depth", "lca_depth_norm", "lca_normalized_height")
      
      #==2 annotation cell state base expression or Fluorescence intensity
      terminal.cell <- data %>% select("cell","time","blot") %>% group_by(cell) %>% 
        mutate(start.time=min(time),end.time=max(time)) %>% ungroup() %>% 
        #time==max(time)
        filter(time==max(time)) %>% select(-time) %>% 
        group_by(cell) %>% mutate(mean.blot.raw=mean(blot)) %>% ungroup() %>% 
        mutate(mean.blot=mean.blot.raw-min(mean.blot.raw)) %>% select(-blot) %>% unique() %>%
        mutate(state=case_when(mean.blot>(max(mean.blot)*0.3)~"ON",mean.blot<=(max(mean.blot)*0.3)~"OFF"))
      #table(terminal.cell$state)
      terminal.cell <- terminal.cell %>% mutate(cell_ID=cell, cellNum=1) %>% select(cell, cell_ID, state, cellNum)
      colnames(terminal.cell) <- c("nodeLabel", "cell_ID", "celltype", "cellNum")
      ##
      cell_type <- c("OFF","ON")
      Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
      Bound_matrix[upper.tri(Bound_matrix)] <- 0
      
      ##get simulation
      df_random <- lapply(seq(n), function(i){
        #i <- 1
        repeat{
          cat("times", i, "\n")
          ##random tree leaf
          sim.terminal.cell <- terminal.cell %>% mutate(celltype=sample(.$celltype, size=nrow(.), replace=F))
          #table(terminal.cell$celltype)
          df_depth_try <- tryCatch({
            lapply(seq(length(depth)), function(d){
              #d <- 1
              sel.depth <- depth[d]
              print(sel.depth)
              #== calculate_transMatrix
              shuffle_optimal_results_d <- get_optimal_transition_matrix(
                node_pair_depth = treeDist,
                cell_info = sim.terminal.cell,
                Sel_u = use.depth,
                fi_depth = sel.depth,
                Bound_Matrix = Bound_matrix,
                cell_type_list = cell_type,
                mc.cores = 80
              )
              shuffle_optimal_results_d <- shuffle_optimal_results_d$optimal_norm_df_dataframe %>% 
                as.data.frame() %>% filter(type_combined_new == "OFF_ON")
              
              return(shuffle_optimal_results_d)
            }) %>% bind_rows() %>% mutate(times=i)
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
      df_random <- df_random %>% as.data.frame() %>% mutate(gene=sel.gene, file=sel.file)
      
      ##
      out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/shuffle_allDep_perGene1000/"
      out_filename <- gsub(".csv","",sel.file)
      saveRDS(df_random, file= paste0(out_path, out_filename, ".shuffle", ".Rds"))
      print("Successful")
      
      return(df_random)
    }) %>% bind_rows()
    return(sim.op.trans.sel.gene)
  }) %>% bind_rows()
  return(sim.op.trans.pergene)
}


#============Run===============================================================================
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/filter_op_vs_sheer_branch/Celegans_filter_fileinfo.Rda")
#== 1 get transition rate accuracy base on actual phylo
#== 1.1 get optimal transition by PASTRI
#calculate_lca_depths
all.treeDist <- treeDist.pergene(file_infos=filter_file)
save(all.treeDist,file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/Celegans_treeDist_file_info.Rda")

#load some files
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/Celegans_treeDist_file_info.Rda")
#get optimal trans
mrca.depth.norm_optimal_res_d0.2 <- get.optimal.irreversible.transT.allgene(file_infos = filter_file,
                                                                            all_treeDist = all.treeDist,
                                                                            use.depth = "lca_depth_norm",
                                                                            sel.depth = 0.2) %>% mutate(condition="d0.2")
mrca.depth.norm_optimal_res_d0.4 <- get.optimal.irreversible.transT.allgene(file_infos = filter_file,
                                                                            all_treeDist = all.treeDist,
                                                                            use.depth = "lca_depth_norm",
                                                                            sel.depth = 0.4) %>% mutate(condition="d0.4")
mrca.depth.norm_optimal_res_d0.6 <- get.optimal.irreversible.transT.allgene(file_infos = filter_file,
                                                                            all_treeDist = all.treeDist,
                                                                            use.depth = "lca_depth_norm",
                                                                            sel.depth = 0.6) %>% mutate(condition="d0.6")
mrca.depth.norm_optimal_res_d0.8 <- get.optimal.irreversible.transT.allgene(file_infos = filter_file,
                                                                            all_treeDist = all.treeDist,
                                                                            use.depth = "lca_depth_norm",
                                                                            sel.depth = 0.8) %>% mutate(condition="d0.8")
mrca.depth.norm_optimal_res_d1.0 <- get.optimal.irreversible.transT.allgene(file_infos = filter_file,
                                                                            all_treeDist = all.treeDist,
                                                                            use.depth = "lca_depth_norm",
                                                                            sel.depth = 1.0) %>% mutate(condition="d1.0")
save(mrca.depth.norm_optimal_res_d0.2, mrca.depth.norm_optimal_res_d0.4, mrca.depth.norm_optimal_res_d0.6,
     mrca.depth.norm_optimal_res_d0.8, mrca.depth.norm_optimal_res_d1.0,
     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/Celegans_depth_norm.various.op.irreversible.trans.allgene.Rda")


#== 1.2 get sheer (real) transition rates
depth <- c(0.2, 0.4, 0.6, 0.8, 1.0)
sheer.trans.res <- lapply(seq(length(depth)),function(d){
  #d <- 1
  sel_depth <- depth[d]
  ##set output file path
  sheer.trans.results <- paste0("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/", "sheer.trans.mrca.depth.norm_d", sel_depth, ".Rds")
  ##
  sheer.trans.results_d <- get.sheer.trans.branch(file_infos = filter_file,
                                                  use.depth = "lca_depth_norm",
                                                  sel.depth = sel_depth)
  ## save results
  saveRDS(sheer.trans.results_d, file = sheer.trans.results)
  return(sheer.trans.results_d)
})

#== 1.3 compare PASTRI vs sheer trans results
rm(list=ls())
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/Celegans_depth_norm.various.op.irreversible.trans.allgene.Rda")
depth <- c(0.2, 0.4, 0.6, 0.8, 1.0)

op.vs.sheer.real <- mclapply(seq(length(depth)), function(d){
  #d <- 2
  #sel.depth <- 0.2
  sel.depth <- depth[d]
  #==1 load sheer trans file
  sheer_file_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/"
  sheer_file <- paste0(sheer_file_path,"sheer.trans.mrca.depth.norm_d", sel.depth ,".Rds")
  sheer.trans.results_d <- readRDS(sheer_file)
  sheer.trans.results_d <- sheer.trans.results_d %>% mutate(name=paste(type.combined, gene, file, sep = "_")) %>%
    select(name, sheer.T)
  #==2 get optimal trans
  sel.op.trans <- paste0("mrca.depth.norm_optimal_res_d", sprintf("%.1f", sel.depth)) %>% get()
  
  #==3 merge optimal and sheer trans probability
  op_d <- sel.op.trans %>% select(type_combined_new, gene, file, rep, norm_optimal_T) %>%
    mutate(name=paste(type_combined_new, gene, file, sep = "_")) %>% unique() %>%
    mutate(norm_optimal_T = ifelse(norm_optimal_T < 0, 0, norm_optimal_T))
  op.sheer_d <- full_join(op_d, sheer.trans.results_d, by="name") %>% filter(!is.na(norm_optimal_T)) %>%
    select(type_combined_new, gene, file, rep, norm_optimal_T, sheer.T) %>%
    mutate(sheer.T = ifelse(is.na(sheer.T), 0, sheer.T)) %>% 
    filter(type_combined_new=="OFF_ON")
  
  #==4 get mean trans
  gene <- unique(op.sheer_d$gene)
  op.sheer.mean <- mclapply(seq(length(gene)),function(ge){
    #ge <- 1
    sel.gene <- gene[ge]
    trans.mean <- op.sheer_d %>% filter(gene==sel.gene) %>% group_by(type_combined_new) %>%
      mutate(op.mean=mean(norm_optimal_T), op.sd=sd(norm_optimal_T), op.se=op.sd/sqrt(n()),
             sheer.mean=mean(sheer.T), sheer.sd=sd(sheer.T), sheer.se=sheer.sd/sqrt(n())) %>%
      mutate(op.sd=ifelse(is.na(op.sd), 0, op.sd), op.se=ifelse(is.na(op.se), 0, op.se),
             sheer.sd=ifelse(is.na(sheer.sd), 0, sheer.sd), sheer.se=ifelse(is.na(sheer.se), 0, sheer.se)) %>%
      as.data.frame()
    return(trans.mean)
  },mc.cores = 60) %>% bind_rows() %>% select(-c(norm_optimal_T, sheer.T,file,rep)) %>% distinct()
  
  return(data.frame(cor=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "pearson")$estimate,
                    cor.p=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "pearson")$p.value,
                    rho=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "spearman")$estimate,
                    rho.p=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "spearman")$p.value,
                    fdist=dist(rbind(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean), method = "euclidean") %>% as.vector(),
                    mhdist=dist(rbind(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean), method = "manhattan") %>% as.vector(),
                    times="real",
                    depth=sel.depth,
                    stringsAsFactors =F))
}) %>% bind_rows()

#op.vs.sheer.real
saveRDS(op.vs.sheer.real, file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/op.vs.sheer.real.Rds")

##
plot.op.vs.sheer <- lapply(seq(length(depth)), function(d){
  #d <- 2
  #sel.depth <- 0.2
  sel.depth <- depth[d]
  print(sel.depth)
  #==1 load sheer trans file
  sheer_file_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/"
  sheer_file <- paste0(sheer_file_path,"sheer.trans.mrca.depth.norm_d", sel.depth ,".Rds")
  sheer.trans.results_d <- readRDS(sheer_file)
  sheer.trans.results_d <- sheer.trans.results_d %>% mutate(name=paste(type.combined, gene, file, sep = "_")) %>%
    select(name, sheer.T)
  #==2 get optimal trans
  sel.op.trans <- paste0("mrca.depth.norm_optimal_res_d", sprintf("%.1f", sel.depth)) %>% get()
  
  #==3 merge optimal and sheer trans probability
  op_d <- sel.op.trans %>% select(type_combined_new, gene, file ,rep, norm_optimal_T) %>%
    mutate(name=paste(type_combined_new, gene, file, sep = "_")) %>% unique() %>%
    mutate(norm_optimal_T = ifelse(norm_optimal_T < 0, 0, norm_optimal_T))
  op.sheer_d <- full_join(op_d, sheer.trans.results_d, by="name") %>% filter(!is.na(norm_optimal_T)) %>%
    select(type_combined_new, gene, file, rep, norm_optimal_T, sheer.T) %>%
    mutate(sheer.T = ifelse(is.na(sheer.T), 0, sheer.T)) %>% 
    filter(type_combined_new=="OFF_ON")
  
  #==4 get mean trans
  gene <- unique(op.sheer_d$gene)
  op.sheer.mean <- mclapply(seq(length(gene)),function(ge){
    #ge <- 1
    sel.gene <- gene[ge]
    trans.mean <- op.sheer_d %>% filter(gene==sel.gene) %>% group_by(type_combined_new) %>%
      mutate(op.mean=mean(norm_optimal_T), op.sd=sd(norm_optimal_T), op.se=op.sd/sqrt(n()),
             sheer.mean=mean(sheer.T), sheer.sd=sd(sheer.T), sheer.se=sheer.sd/sqrt(n())) %>%
      mutate(op.sd=ifelse(is.na(op.sd), 0, op.sd), op.se=ifelse(is.na(op.se), 0, op.se),
             sheer.sd=ifelse(is.na(sheer.sd), 0, sheer.sd), sheer.se=ifelse(is.na(sheer.se), 0, sheer.se)) %>%
      as.data.frame()
    return(trans.mean)
  },mc.cores = 60) %>% bind_rows() %>% select(-c(norm_optimal_T, sheer.T,file,rep)) %>% distinct()
  
  res <- data.frame(cor=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "pearson")$estimate,
                    cor.p=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "pearson")$p.value,
                    rho=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "spearman")$estimate,
                    rho.p=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "spearman")$p.value,
                    fdist=dist(rbind(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean), method = "euclidean") %>% as.vector(),
                    mhdist=dist(rbind(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean), method = "manhattan") %>% as.vector(),
                    times="real",
                    depth=sel.depth,
                    stringsAsFactors =F)
  #==5 plot
  highlighted_point <- op.sheer.mean[op.sheer.mean$gene == "mep-1", ]
  #"CD20080802_die-1_sd1566.csv" "CD20080814_lin-13_6B12_1_L2.csv" "CD20090401_mep-1_11C1_8_L1.csv"
  if(sel.depth==0.4){
    p <- ggplot(op.sheer.mean, aes(x=sheer.mean, y=op.mean))+geom_point(size=3,color="#FE9999")+theme_classic()+
      labs(title=paste0("d_",sel.depth,
                        "\ncor = ", round(res$cor, 2), " P=", round(log10(res$cor.p),3),
                        #"\nrho = ", round(res$rho, 3), " P=", res$rho.p,
                        "\neuclidean distance = ", round(res$fdist, 2)),
           #"\nmanhattan distance = ", round(res$mhdist, 2)),
           x="Average sheer switching probability",y="Average inferred switching probability")+
      geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+#scale_shape_manual(values=c(17,15,16,18))+
      geom_errorbar(aes(xmin=sheer.mean-sheer.se,xmax=sheer.mean+sheer.se),color="gray",width=0)+xlim(0,1)+
      geom_errorbar(aes(ymin=op.mean-op.se,ymax=op.mean+op.se),color="gray",width=0)+
      geom_point(data=highlighted_point, aes(x=sheer.mean, y=op.mean), color="#1F78B4", size=5, shape=19)+
      theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
            legend.title=element_blank(),legend.text=element_text(size=10,color="black"))
  }else{
    p <- ggplot(op.sheer.mean, aes(x=sheer.mean, y=op.mean))+geom_point(size=3,color="#FE9999")+theme_classic()+
      labs(title=paste0("d_",sel.depth,
                        "\ncor = ", round(res$cor, 2), " P=", round(log10(res$cor.p),3),
                        #"\nrho = ", round(res$rho, 3), " P=", res$rho.p,
                        "\neuclidean distance = ", round(res$fdist, 2)),
           #"\nmanhattan distance = ", round(res$mhdist, 2)),
           x="Average sheer switching probability",y="Average inferred switching probability")+
      geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+#scale_shape_manual(values=c(17,15,16,18))+
      geom_errorbar(aes(xmin=sheer.mean-sheer.se,xmax=sheer.mean+sheer.se),color="gray",width=0)+xlim(0,1)+
      geom_errorbar(aes(ymin=op.mean-op.se,ymax=op.mean+op.se),color="gray",width=0)+
      #geom_point(data=highlighted_point, aes(x=op.mean, y=sheer.mean), color="#1F78B4", size=5, shape=19)+
      theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
            legend.title=element_blank(),legend.text=element_text(size=10,color="black"))
  }
  
  return(p)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure3/Fig3C.pdf",width=4.6,height=5)  
plot.op.vs.sheer[[2]]
dev.off()


#=====================================================================================================
#== 2 get transition rate accuracy base on shuffle phylo
#== 2.1 PASTRI-inferred trans of shuffle phylo
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/filter_op_vs_sheer_branch/Celegans_filter_fileinfo.Rda")
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/Celegans_treeDist_file_info.Rda")

df_radom_allgene <- get.sim.op.irreversible.transT.allgene(file_infos = filter_file,
                                                           all_treeDist = all.treeDist,
                                                           use.depth = "lca_depth_norm",
                                                           depth = c(0.2, 0.4, 0.6, 0.8, 1.0),
                                                           n=1000)

out_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/"
saveRDS(df_radom_allgene, file= paste0(out_path, "1000shuffle_allDep",".Rds")) #c(0.2, 0.4, 0.6, 0.8, 1.0)

#== 2.2 compare sheer trans of real phylo and PASTRI-inferred trans of shuffle phylo
rm(list=ls())
path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/"
sim.op.trans <- readRDS(paste0(path, "1000shuffle_allDep",".Rds"))

depth <- c(0.2, 0.4, 0.6, 0.8, 1.0)
n <- 1000

df_random <- lapply(seq(n), function(i){
  #i <- 1
  cat("times", i, "\n")
  sel.sim.op.trans_1 <- sim.op.trans %>% filter(times==i)
  
  df_random_depth <- lapply(seq(length(depth)), function(d){
    #d <- 2
    sel.depth <- depth[d]
    print(sel.depth)
    #==sim
    sel.sim.op.trans_2 <- sel.sim.op.trans_1 %>% filter(depth==paste0("d_", sel.depth)) %>%
      mutate(norm_optimal_T = ifelse(norm_optimal_T < 0, 0, norm_optimal_T))
    #==sheer
    sheer_file_path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/"
    sheer_file <- paste0(sheer_file_path, "sheer.trans.mrca.depth.norm_d", sel.depth ,".Rds")
    sheer.trans.results_d <- readRDS(sheer_file) %>% filter(type.combined=="OFF_ON") %>%
      select(file, sheer.T)
    
    #==merge sim optimal and sheer trans probability
    sel.sim.sheer <- full_join(sel.sim.op.trans_2, sheer.trans.results_d, by="file") %>%
      select(depth, times, gene, file, type_combined_new, norm_optimal_T, sheer.T) %>%
      filter(!is.na(norm_optimal_T)) %>% mutate(sheer.T = ifelse(is.na(sheer.T), 0, sheer.T))
    #cor.test(sel.sim.sheer$norm.optimal.T, sel.sim.sheer$sheer.T, method="pearson")
    
    #==get mean trans
    gene <- unique(sel.sim.sheer$gene)
    sel.sim.sheer.mean <- lapply(seq(length(gene)), function(ge){
      #ge <- 1
      sel.gene <- gene[ge]
      trans.mean <- sel.sim.sheer %>% filter(gene==sel.gene) %>% group_by(type_combined_new) %>%
        mutate(op.mean=mean(norm_optimal_T), sheer.mean=mean(sheer.T)) %>% ungroup() %>% as.data.frame()
      return(trans.mean)
    }) %>% bind_rows() %>%
      select(depth, times, gene, type_combined_new, op.mean, sheer.mean) %>% unique()
    
    #==output
    res.cor <- cor.test(sel.sim.sheer.mean$op.mean, sel.sim.sheer.mean$sheer.mean, method = "pearson")
    res.rho <- cor.test(sel.sim.sheer.mean$op.mean, sel.sim.sheer.mean$sheer.mean, method = "spearman")
    res.fdist <- dist(rbind(sel.sim.sheer.mean$op.mean, sel.sim.sheer.mean$sheer.mean), method = "euclidean") %>% as.vector()
    res.mhdist <- dist(rbind(sel.sim.sheer.mean$op.mean, sel.sim.sheer.mean$sheer.mean), method = "manhattan") %>% as.vector()
    
    return(data.frame(depth=sel.depth,
                      times=i,
                      cor=res.cor$estimate, cor.pval=res.cor$p.value,
                      rho=res.rho$estimate, rho.pval=res.rho$p.value,
                      fdist=res.fdist, mhdist=res.mhdist,
                      stringsAsFactors =F))
  }) %>% bind_rows()
  return(df_random_depth)
}) %>% bind_rows()

saveRDS(df_random, file= paste0(path, "1000allDep_compare_shuffle_allLeafs_sim.vs.sheer",".Rds")) #c(0.2, 0.4, 0.6, 0.8, 1.0)

#== 2.3 get PASTRI-inferred accuracy
rm(list=ls())
path <- "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/"
real.res <- readRDS(paste0(path, "op.vs.sheer.real.Rds"))
sim.res <- readRDS(paste0(path, "1000allDep_compare_shuffle_allLeafs_sim.vs.sheer.Rds"))

depth <- c(0.2, 0.4, 0.6, 0.8, 1.0)

#get p value base on rank
sim.real.res_p_rank <- lapply(seq(length(depth)), function(d){
  #d <- 1
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>% 
    select(depth, times, cor, fdist, mhdist)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(which(!grepl("pval", names(.)))) %>%
    select(depth, times, cor, fdist, mhdist)
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
  filter(methods=="fdist") %>% select(-methods)
colnames(sim.real.res_p_rank) <- c("depth", "p_rank")

#get p value base on test
sim.real.res_p.new <- lapply(seq(length(depth)), function(d){
  #d <- 2
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(depth,times,cor,rho,fdist,mhdist) #%>%filter(times>100&times<=201)
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(depth,times,cor,rho,fdist,mhdist)
  ##
  diff <- sel.op_real$fdist-mean(sel.op_shuffle$fdist)
  z_score <- (sel.op_real$fdist - mean(sel.op_shuffle$fdist)) / sd(sel.op_shuffle$fdist)
  ##
  #cor.p <- wilcox.test(sel.op_shuffle$cor, mu=sel.op_real$cor)$p.value %>% -log(.)
  #fdist.p <- wilcox.test(sel.op_shuffle$fdist, mu=sel.op_real$fdist)$p.value %>% -log(.)
  #mhdist.p <- wilcox.test(sel.op_shuffle$mhdist, mu=sel.op_real$mhdist)$p.value %>% -log10(.)
  #options(digits = 22)
  #cor.p <- t.test(sel.op_shuffle$cor, mu=sel.op_real$cor)$p.value %>% -log(.)
  #fdist.p <- t.test(sel.op_shuffle$fdist, mu=sel.op_real$fdist)$p.value %>% -log(.)
  fdist.t.p <- t.test(sel.op_shuffle$fdist, mu=sel.op_real$fdist)$statistic %>%
    pt(., df=999, lower.tail=F, log.p=T)
  fdist.t.p <- -((fdist.t.p/log(10))+log10(2))
  #wilcox.test
  fdist.wilcox.p <- wilcox.test(sel.op_shuffle$fdist, mu=sel.op_real$fdist,log.p=TRUE)$p.value # %>% -log10(.)
  #a <- t.test(sel.op_shuffle$fdist[1:10], mu=1.9)$statistic %>% pt(., df=9, lower.tail=F, log.p=T)
  #-((a/log(10))+log10(2))
  norm.p.ks <- ks.test(sel.op_shuffle$fdist,"pnorm")$p.value
  norm.p.shapiro <- shapiro.test(sel.op_shuffle$fdist)$p.value
  norm.p.ad <- ad.test(sel.op_shuffle$fdist)$p.value
  
  return(data.frame(depth = sel.depth, Diff = diff, z_score = z_score,
                    #real.cor=sel.op_real$cor, cor.Pvalue=cor.p,
                    real.fdist=sel.op_real$fdist, fdist.t.Pvalue=fdist.t.p, #t.value=t.value,norm.p.new=norm.p.new,
                    norm.p.ks=norm.p.ks,
                    norm.p.shapiro=norm.p.shapiro,
                    norm.p.ad=norm.p.ad,fdist.wilcox.Pvalue=fdist.wilcox.p,
                    #real.mhdist=sel.op_real$mhdist, mhdist.Pvalue=mhdist.p,
                    stringsAsFactors = F))
}) %>% bind_rows() %>%
  mutate(condition=paste0("<=",depth)) %>% select(-c(Diff, norm.p.ks))


load(paste0(path, "Celegans_depth_norm.various.op.irreversible.trans.allgene.Rda"))
op.vs.sheer.real.bootstrap <- mclapply(seq(length(depth)), function(d){
  #d <- 2
  #sel.depth <- 0.2
  sel.depth <- depth[d]
  #==1 load sheer trans file
  sheer_file <- paste0(path,"sheer.trans.mrca.depth.norm_d", sel.depth ,".Rds")
  sheer.trans.results_d <- readRDS(sheer_file)
  sheer.trans.results_d <- sheer.trans.results_d %>% mutate(name=paste(type.combined, gene, file, sep = "_")) %>%
    select(name, sheer.T)
  #==2 get optimal trans
  sel.op.trans <- paste0("mrca.depth.norm_optimal_res_d", sprintf("%.1f", sel.depth)) %>% get()
  
  #==3 merge optimal and sheer trans probability
  op_d <- sel.op.trans %>% select(type_combined_new, gene, file, rep, norm_optimal_T) %>%
    mutate(name=paste(type_combined_new, gene, file, sep = "_")) %>% unique() %>%
    mutate(norm_optimal_T = ifelse(norm_optimal_T < 0, 0, norm_optimal_T))
  op.sheer_d <- full_join(op_d, sheer.trans.results_d, by="name") %>% filter(!is.na(norm_optimal_T)) %>%
    select(type_combined_new, gene, file, rep, norm_optimal_T, sheer.T) %>%
    mutate(sheer.T = ifelse(is.na(sheer.T), 0, sheer.T)) %>% 
    filter(type_combined_new=="OFF_ON")
  
  #==4 get mean trans
  gene <- unique(op.sheer_d$gene)
  op.sheer.mean <- mclapply(seq(length(gene)),function(ge){
    #ge <- 1
    sel.gene <- gene[ge]
    trans.mean <- op.sheer_d %>% filter(gene==sel.gene) %>% group_by(type_combined_new) %>%
      mutate(op.mean=mean(norm_optimal_T), op.sd=sd(norm_optimal_T), op.se=op.sd/sqrt(n()),
             sheer.mean=mean(sheer.T), sheer.sd=sd(sheer.T), sheer.se=sheer.sd/sqrt(n())) %>%
      mutate(op.sd=ifelse(is.na(op.sd), 0, op.sd), op.se=ifelse(is.na(op.se), 0, op.se),
             sheer.sd=ifelse(is.na(sheer.sd), 0, sheer.sd), sheer.se=ifelse(is.na(sheer.se), 0, sheer.se)) %>%
      as.data.frame()
    return(trans.mean)
  },mc.cores = 60) %>% bind_rows() %>% 
    select(-c(norm_optimal_T, sheer.T,file,rep)) %>% distinct()
  
  res <-data.frame(cor=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "pearson")$estimate,
                   cor.p=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "pearson")$p.value,
                   rho=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "spearman")$estimate,
                   rho.p=cor.test(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean, method = "spearman")$p.value,
                   fdist=dist(rbind(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean), method = "euclidean") %>% as.vector(),
                   mhdist=dist(rbind(op.sheer.mean$op.mean, op.sheer.mean$sheer.mean), method = "manhattan") %>% as.vector(),
                   times="real",
                   depth=sel.depth,
                   stringsAsFactors =F)
  
  #==5 bootstrap sd 
  op.sheer.mean.bootstrap <- lapply(seq(1000), function(i){
    op.sheer.mean.new <- op.sheer.mean[sample(nrow(op.sheer.mean), replace = TRUE), ]
    res1 <- data.frame(cor=cor.test(op.sheer.mean.new$op.mean, op.sheer.mean.new$sheer.mean, method = "pearson")$estimate,
                       rho=cor.test(op.sheer.mean.new$op.mean, op.sheer.mean.new$sheer.mean, method = "spearman")$estimate,
                       fdist=dist(rbind(op.sheer.mean.new$op.mean, op.sheer.mean.new$sheer.mean), method = "euclidean") %>% as.vector(),
                       mhdist=dist(rbind(op.sheer.mean.new$op.mean, op.sheer.mean.new$sheer.mean), method = "manhattan") %>% as.vector(),
                       times=i,
                       depth=sel.depth,
                       stringsAsFactors =F)
    return(res1)
  }) %>% bind_rows()
  
  op.sheer.mean.bootstrap.new <- op.sheer.mean.bootstrap %>%
    mutate(cor.sd=sd(cor), rho.sd=sd(rho), fdist.sd=sd(fdist), mhdist.sd=sd(mhdist)) %>% 
    select(cor.sd, rho.sd, fdist.sd, mhdist.sd, depth) %>% unique()
  
  res2 <- full_join(res, op.sheer.mean.bootstrap.new, by="depth")
  
  return(res2)
}) %>% bind_rows() %>%
  mutate(condition=paste0("<=",depth)) %>% select(condition, fdist, fdist.sd)

sim.real.res_p.new <- full_join(sim.real.res_p.new, op.vs.sheer.real.bootstrap, by="condition")
sim.real.res_p.new <- full_join(sim.real.res_p.new, sim.real.res_p_rank, by="depth")

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure3/Fig3D.pdf",width=4.8,height=4)
ggplot(sim.real.res_p.new, aes(x=condition,y=real.fdist,color=condition))+geom_point(size=3)+
  geom_text(aes(label = p_rank), vjust = -1, size = 3, color = "black") +
  labs(x="Normalized depth (d)",y="Accuracy of PASTRI inference\nEuclidean distance between true and inferred rates")+
  geom_errorbar(aes(ymin=real.fdist-fdist.sd, ymax=real.fdist+fdist.sd),width=0)+
  theme_classic()+#scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values = ifelse(sim.real.res_p.new$condition == "<=0.4", "#91C893", "#D6D6D6"))+
  theme(plot.title = element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"), 
        axis.text=element_text(size=10,color="black"),axis.line=element_line(color="black"), 
        panel.spacing = unit(1, "lines"), legend.position = "none")
dev.off()

##Fig3E
rm(list = setdiff(ls(), c("path", "real.res", "sim.res", "depth", "sim.real.res_p.new")))

plot.shuffle <- lapply(seq(length(depth)), function(d){
  #d <- 2
  sel.depth <- depth[d]
  print(sel.depth)
  #
  sel.op_shuffle <- sim.res %>% filter(depth==sel.depth) %>% select(depth, fdist) 
  #
  sel.op_real <- real.res %>% filter(depth==sel.depth) %>% select(depth, fdist)
  #
  sel.sim.real.res_p <- sim.real.res_p.new %>% filter(depth==sel.depth)
  
  cutoff_1 <- sel.op_real$fdist+0.1 
  cutoff_2 <- min(sel.op_shuffle$fdist)
  
  #plot
  p <- ggplot(sel.op_shuffle, aes(x=fdist))+ geom_histogram(binwidth = 0.0002,fill="#92C5DE")+
    #geom_freqpoly(binwidth=0.001,linewidth=0.1)
    #geom_freqpoly(bins=2000,col="#1D91C0")+
    geom_vline(data = sel.op_real, aes(xintercept = fdist), linetype = "dashed", color="#D73027")+
    labs(title="Permutation test P = 0.001", x="", y="Count")+
    xlim(cutoff_1-0.2, max(sel.op_shuffle$fdist)+0.001)+
    scale_x_break(c(cutoff_1,cutoff_2),scales = 5,space=0.2)+theme_classic()+
    scale_y_continuous(expand=c(0,0))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),
          panel.spacing = unit(1, "lines"),plot.title=element_text(size=10,color="black",hjust=0.5),
          axis.title=element_text(size=10,color="black"),axis.text=element_text(size=10,color="black"),
          axis.line=element_line(color="black"),axis.line.x.top=element_blank(),
          axis.text.x.top = element_blank(),axis.ticks.x.top=element_blank())
  return(p)
})

pdf("/mnt/data5/disk/yangwj/Result_plots/Figure3/Fig3E.pdf",width=5,height=4.5)  
plot.shuffle[[2]]
dev.off()

