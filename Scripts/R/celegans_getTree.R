rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  library(parallel)
  library(dplyr)
  library(igraph)
  library(ape)
  library(ggtree)
  library(ggplot2)
})

setwd("/mnt/data5/disk/yangwj/Scripts/Cell_State_Transition/Celegans_EPIC")
## ==================== source and define functions ===========================================
trans_linstr_to_phylo <- function(linstr){
  ## 1. pre-define specific string relationship
  specific_parent_infos <-
    data.frame(from = c("P0","P0","P1", "P1", "EMS", "EMS", "P2", "P2", "P3", "P3", "P4", "P4","AB","AB"),
               to = c("AB","P1","EMS", "P2", "MS", "E", "C", "P3", "D", "P4", "Z2", "Z3","ABa","ABp"),
               stringsAsFactors =F)
  ## 2. contruct the graph df
  #linstr <- all_nodes
  net_df <- 
    data.frame(from = substr(linstr, start = 0, stop = nchar(linstr) - 1),
               to = linstr,
               stringsAsFactors = F)
  #net_df <- filter(net_df, from %in% linstr)
  #net_df %>% filter(!from %in% linstr)
  from.node.fi <- c("A","P","EM","M","","Z","Nuc")
  net_df <- net_df %>% filter(!from %in% from.node.fi)
  net_df <- bind_rows(specific_parent_infos, net_df) %>% unique()
  ## 3. from graph to phylo
  g <- graph_from_data_frame(net_df, directed = TRUE)
  #edges <- get.edgelist(g)
  tree_obj <- as.phylo(g)
  tree_obj$edge.length <- rep(1, nrow(tree_obj$edge)) 
  return(tree_obj)
}

tree.pergene <- function(file_infos){
  gene_list <- unique(file_infos$gene)
  res.tree.pergene <- mclapply(seq(length(gene_list)),function(g){
    #g <-1
    sel.gene <- gene_list[g]
    sel.fileInfo <- file_infos %>% filter(gene==sel.gene)
    file_list_pergene <- unique(sel.fileInfo$filename)
    res.sel.gene <- mclapply(seq(length(file_list_pergene)),function(f){
      #f <-1
      sel.file <- file_list_pergene[f]
      data <- read.csv(sel.file,stringsAsFactors=FALSE)
      all_nodes <- unique(data$cell) %>% unlist()
      trans_tree <- trans_linstr_to_phylo(all_nodes)
      #output_filename <- paste0(sel.gene, "_", basename(sel.file), ".nwk")
      output_filename <- paste0(basename(sel.file), ".nwk")
      output_path <- file.path("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/treefile_Celegans", output_filename)
      write.tree(trans_tree, file = output_path)
      return(output_path)
      #return(trans_tree)
    },mc.cores = 80)
    return(unlist(res.sel.gene))
  })
  return(unlist(res.tree.pergene))
}

get.numLeaf.pergene <- function(file_infos){
  #file_infos=file_info
  gene_list <- unique(file_infos$gene)
  res.tree.pergene <- mclapply(seq(length(gene_list)),function(g){
    #g <- 1
    sel.gene <- gene_list[g]
    sel.fileInfo <- file_infos %>% filter(gene==sel.gene)
    file_list_pergene <- unique(sel.fileInfo$filename)
    res.sel.gene <- mclapply(seq(length(file_list_pergene)),function(f){
      #f <-1
      sel.file <- file_list_pergene[f]
      data <- read.csv(sel.file,stringsAsFactors=FALSE)
      #== get endpoint cell number
      terminal.cell <- data  %>% select(cell,time) %>% filter(time==max(time))
      #
      res <- data.frame(gene=sel.gene,
                        filename=sel.file,
                        num_terminal_cell=length(unique(terminal.cell$cell)),
                        stringsAsFactors=FALSE)
      return(res)
    }, mc.cores = 60) %>% bind_rows()
    return(res.sel.gene)
  }) %>% bind_rows()
  return(res.tree.pergene)
}


#============Run===============================================================================
file_info <- read.csv("/mnt/data5/disk/yangwj/Scripts/Cell_State_Transition/Celegans_EPIC/A_fileInfos.new.csv",stringsAsFactors=FALSE)
file_info <- file_info %>% group_by(gene) %>% add_count(gene,name="rep") %>% ungroup()

#====== 0 get phylo tree file (nwk) ======
tree_pergene <- tree.pergene(file_infos = file_info)

#====== 1 get terminal cell nums per tree base on csv file ======
numLeaf.pergene <- get.numLeaf.pergene(file_infos=file_info)
filter_file <- numLeaf.pergene %>% filter(num_terminal_cell > 200) %>% 
  group_by(gene) %>% add_count(gene, name="rep") %>% ungroup()
save(filter_file, file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig3_Celegans/filter_op_vs_sheer_branch/Celegans_filter_fileinfo.Rda")


