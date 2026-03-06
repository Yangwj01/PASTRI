rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  library(tidyverse)
  library(parallel)
  library(Seurat)
  library(reshape2)
  library(magrittr)
  library(ggtree)
  library(ggplot2)
  library(treeio)
  library(ape)
  library(patchwork)
  library(RColorBrewer)
})

## ==================== source and define functions ===========================================
## depth: distance to root 
get.node.depths <- function(info.df){
  root <- info.df %>% filter(., from == "0" | to == from) %>% `[`("from") %>% unlist() %>% unique()
  node.depth.df <- data.frame(to = root, depth = 0, stringsAsFactors = F)
  current.from <- root
  d <- 1
  repeat{
    source.df <- filter(info.df, from %in% current.from)
    if (nrow(source.df) == 0) {
      break
    } 
    current.from <- unique(source.df$to)
    re.df <- data.frame(to = as.character(source.df$to %>% unlist()),depth = d,stringsAsFactors = F)
    node.depth.df <- bind_rows(node.depth.df, re.df)
    d <- d + 1
  }
  node.depth.df$to <- as.character(node.depth.df$to)
  return(node.depth.df)
}

## distance to bottom
get.node.hieghts <- function(info.df){
  all.inodes <- unique(info.df$from)
  all.inodes <- all.inodes[all.inodes != "0"]
  
  node.depth.df <- 
    mclapply(all.inodes, mc.cores=40, function(n){
      d <- 1
      use.node <- n
      repeat{
        source.df <- filter(info.df, from %in% use.node)
        if (nrow(source.df) == 0) {
          break
        } 
        d <- d + 1
        use.node <- source.df$to %>% unlist() 
      }
      return(data.frame(to = n,
                        height = d))
    }) %>% bind_rows()
  node.depth.df$to <- as.character(node.depth.df$to)
  return(node.depth.df)
}

sample.tree.infos <- function(bc.cell, sample.name, tree.path, alleleinfo.path){
  rawtree <- read.tree(tree.path)
  alleleinfo <- read.csv(alleleinfo.path, stringsAsFactors = F)
  alleleinfo$BC <- gsub("BC_","",alleleinfo$BC)
  #colnames(alleleinfo) <- c("node","BC","editMutation","cellNum","sample","celltype")
  # extract sample bc 
  sub.bc.ct <- bc.cell %>% filter(., sample == sample.name)
  # map celltype to tree node
  allele.ct <- left_join(alleleinfo, sub.bc.ct) %>% 
    dplyr::select(node, BC, sample, cellNum,celltype)  
  # join with tree df
  rawtree.df <- as_tibble(rawtree)
  rawtree.info.df <- left_join(rawtree.df, allele.ct) %>% 
    mutate(to = ifelse(label %in% c("", "1"), paste0(sample.name, "_", node), paste0(sample.name, "_", label, "_", BC)),
           from = paste0(sample.name, "_", parent),
           type = ifelse(label %in% c("", "1"), ifelse(label == "", "root", "inode"), "leaf")) %>% 
    dplyr::select(-c(parent, node, branch.length))
  
  rawtree.info.df$sample[is.na(rawtree.info.df$sample)] <- sample.name
  rawtree.info.df$celltype <- as.character(rawtree.info.df$celltype)
  rawtree.info.df$celltype[is.na(rawtree.info.df$celltype)] <- "inter"
  #########################################
  ## =========== optional: modify the tree topology ============================
  ######
  #### ==== 1. for every multi cell node ,split it by cell number and sign them to an additional internal node ==== ###
  message("modify multifurcat nodes")
  new.rawtree.info <- 
    rawtree.info.df %>% group_by(., label) %>% do({
      mydf <- .
      if (mydf$type[1] %in% c("root", "inode")) {
        re.df <- mydf
      } else if (nrow(mydf) == 1) {
        re.df <- mydf
      } else {
        # record sample label and original from 
        sample.label <- mydf$sample[1]
        ori.from <- mydf$from[1]
        # create an additional internal node 
        add.inode.label <- paste0("add_", mydf$label[1])
        add.inode.to <- paste0(sample.label, "_", add.inode.label)
        # add the node to raw df and return
        re.df <- 
          mydf %>% mutate(from = add.inode.to) %>%
          add_row(label = add.inode.label, 
                  sample=sample.label, 
                  celltype="inter",
                  to = add.inode.to,
                  from = ori.from, 
                  type = "inode")
      }
      re.df
    })
  
  ### ==== 2. for every internal node, if it have more than one leafs and have internal node descendants simultaneously ==
  # add an additional internal node as the ancestor of isolated leafs 
  # test
  message("modify internal node")
  new.tree.info <- new.rawtree.info %>% group_by(from) %>% do({
    sub.spl <- .
    if (all(sub.spl$type == "leaf") | (sum(sub.spl$type == "leaf") %in% c(0,1))) {
      new.sub.spl <- sub.spl
    } else {
      # record sample label and original from 
      sample.label <- sub.spl$sample[1]
      ori.from <- sub.spl$from[1]
      # create an additional internal node 
      add.inode.label <- paste0(ori.from, "_add")
      add.inode.to <- paste0(sample.label, "_", add.inode.label)
      # add the node to raw df and return
      new.sub.spl <- 
        sub.spl %>% group_by(type) %>% do({
          mydf <- .
          if (mydf$type[1] %in% c("inode", "root")) {
            re.df <- mydf
          } else {
            re.df <- mydf %>% 
              mutate(from = add.inode.to) %>% 
              add_row(label = add.inode.label, 
                      sample=sample.label, 
                      celltype="inter",
                      to = add.inode.to,
                      from = ori.from, 
                      type = "inode")
          }
          re.df
        })
    }
    new.sub.spl
  })
  
  return(new.tree.info)
}

plot.TreeAndCelltype <- function(tree, all.info.df, sample.name){
  # 1. extract all node celltype infos of all nodes
  node.celltype.info <- all.info.df %>% filter(sample == sample.name & celltype != "inter")
  # 2. note all multi cell nodes
  select.node.celltype <- node.celltype.info %>% ungroup() %>%
    mutate(celltype=ifelse(cellNum>1, "multi_nodes", celltype)) %>%
    select(label, sample, cellNum, celltype, type) %>% unique() %>% as.data.frame()
  rownames(select.node.celltype) <- select.node.celltype$label
  
  # 4. set use colors
  #
  celltype.colors <- structure(c("#B2DF8A","#3DA42D","#FE9999","#A4CEDE","#FDBE6F","#FF9116","#BBBFE0","#673A95","#B25A2C","gray"),
                               names=as.character(c(paste0("C",seq(0,8)), "multi_nodes")))
  # 5. only keep the celltype cols
  select.node.celltype.use <- as.data.frame(select.node.celltype[, "celltype", drop=F])
  #rownames(select.node.celltype.use) <- select.node.celltype$label
  select.node.celltype.use$celltype <- as.factor(select.node.celltype.use$celltype)
  
  # plot circular tree and celltypes
  p <- ggtree(tree,branch.length="none",layout="circular",size=0.1,color="black",alpha=1)+xlim(-10, NA) #circular or radial
  p <- rotate_tree(p, 80)
  p <- 
    gheatmap(p, select.node.celltype.use, offset=0.05, width=0.1, colnames_angle=120, colnames_offset_y = .25, colnames = FALSE,color=NA) + xlim(-2, NA) +
    scale_fill_manual(name="Cell types",values = celltype.colors)
  p <- p+theme(legend.position = "right")

    return(p)
}

## ============Run=============================================================================
## ======= read seurat object ====== 
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")
sce.query.hESC <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.query.scobj.SOT.V4.Rds")

# extract cell barcode(BC), nodename, celltype
bc.cell <- rbind(sce@meta.data[,c(1:11, 18:19)], sce.query.hESC@meta.data) %>% as.data.frame() %>%
  mutate(sample = gsub("10X_", "", orig.ident),sample = gsub("H9_", "hESC_", sample),
         BC = gsub(".+_|-1", "", rownames(.))) %>%
  dplyr::select(-orig.ident)

## ====== modify all sample tree structure ======
## set sample paths
pacbio_path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree"
sample.list <- c("HSC_B7","HSC_E12","HSC_F12","hESC_A7","hESC_C9","hESC_C12","hESC_D3") %>% as.character()

for (sample_name in sample.list) {
  #sample_name <- sample.list[4]
  cut_str <- paste0(sample_name, "_B")
  print(cut_str)
  #assign.name <- paste0(gsub("10X_", "", sample_name), ".tree.info")
  assign.name <- paste0(sample_name, ".tree.info")
  # run
  allele.info_path = file.path(pacbio_path, paste0(cut_str, "/", sample_name, ".AllelesInfo.csv"))
  tree_path = file.path(pacbio_path, paste0(cut_str, "/", sample_name, ".nwk"))
  # run
  sub.tree.info <- 
    sample.tree.infos(bc.cell = bc.cell,
                      sample.name = sample_name,
                      tree.path = tree_path, 
                      alleleinfo.path = allele.info_path)
  assign(assign.name,sub.tree.info)
}


## merge them to one dataframe
allSamples <- c("HSC_B7","HSC_E12","HSC_F12","hESC_A7","hESC_C9","hESC_C12","hESC_D3")
new.tree.allSample <- data.frame(NULL)

for(thisSample in allSamples) {
  new.tree.allSample <- bind_rows(get(paste0(thisSample, ".tree.info")), new.tree.allSample)
}
# add a root name '0' to all tree
new.tree.allSample$from[new.tree.allSample$type == "root"] <- "0"

# add height info
all.tree.height <- get.node.hieghts(new.tree.allSample)
new.tree.allSample$height <- all.tree.height$height[match(new.tree.allSample$to, all.tree.height$to)]
new.tree.allSample$from.height <- all.tree.height$height[match(new.tree.allSample$from, all.tree.height$to)]
new.tree.allSample$from.height[new.tree.allSample$from == "0"] <- max(new.tree.allSample$height, na.rm = T) + 1
# add depth info
all.tree.depth <- get.node.depths(new.tree.allSample)
new.tree.allSample$depth <- all.tree.depth$depth[match(new.tree.allSample$to, all.tree.depth$to)]
new.tree.allSample$from.depth <- all.tree.depth$depth[match(new.tree.allSample$from, all.tree.depth$to)]
new.tree.allSample$from.depth[new.tree.allSample$from == "0"] <- min(new.tree.allSample$depth, na.rm = T) - 1

# ====== extract leaf belongs =======
new.dfTip.long <- new.tree.allSample %>%
  filter(type == "leaf") %>% ## only scan the terminal nodes
  group_by(to) %>%
  do({
    myDf <- .;
    curRow <- myDf;
    cumNode <- c(myDf$from);
    while(curRow$type != "root") {
      cumNode <- c(cumNode, curRow$from)
      curRow <- new.tree.allSample %>% 
        filter(to == curRow$from);
    }
    data.frame(intern = cumNode,tip=myDf$label);
  })
new.dfTip.long <- new.dfTip.long[!duplicated(new.dfTip.long),]
new.dfTip.long$celltype <- new.tree.allSample$celltype[match(new.dfTip.long$to, new.tree.allSample$to)]
any(is.na(new.dfTip.long))

### ====== save results ======
saveRDS(new.tree.allSample, 
        file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/all_tree_dataframe_modify_new.Rds")
saveRDS(new.dfTip.long,
        file = "/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/all_tree_modify_tip_long_new.Rds")


## ======= Figure 4G and Figure S6D plot tree ====== 
#1 get tree and note celltype info
tree.infos <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/all_tree_dataframe_modify_new.Rds")
path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/"
B7.tree <- read.tree(paste0(path, "HSC_B7_B/HSC_B7.treefile"))
E12.tree <- read.tree(paste0(path, "HSC_E12_B/HSC_E12.treefile"))
F12.tree <- read.tree(paste0(path, "HSC_F12_B/HSC_F12.treefile"))

A7.tree <- read.tree(paste0(path, "hESC_A7_B/hESC_A7.treefile"))
C9.tree <- read.tree(paste0(path, "hESC_C9_B/hESC_C9.treefile"))
C12.tree <- read.tree(paste0(path, "hESC_C12_B/hESC_C12.treefile"))
D3.tree <- read.tree(paste0(path, "hESC_D3_B/hESC_D3.treefile"))

p.all.B7 <- plot.TreeAndCelltype(tree=B7.tree,all.info.df=tree.infos,sample.name="HSC_B7")
p.all.E12 <- plot.TreeAndCelltype(tree=E12.tree,all.info.df=tree.infos,sample.name="HSC_E12")
p.all.F12 <- plot.TreeAndCelltype(tree=F12.tree,all.info.df=tree.infos,sample.name="HSC_F12")
p.all.A7 <- plot.TreeAndCelltype(tree=A7.tree,all.info.df=tree.infos,sample.name="hESC_A7")
p.all.C9 <- plot.TreeAndCelltype(tree=C9.tree,all.info.df=tree.infos,sample.name="hESC_C9")
p.all.C12 <- plot.TreeAndCelltype(tree=C12.tree,all.info.df=tree.infos,sample.name="hESC_C12")
p.all.D3 <- plot.TreeAndCelltype(tree=D3.tree,all.info.df=tree.infos,sample.name="hESC_D3")

(p.all.B7|p.all.E12|p.all.F12|plot_spacer())/(p.all.A7|p.all.C9|p.all.C12|p.all.D3)

#2 bootstrap
sample <- c("HSC_B7","HSC_E12","HSC_F12", "hESC_A7","hESC_C9","hESC_C12","hESC_D3")
p.bootstrap <- lapply(seq(sample), function(s){
  #s <- 1
  sel.sample <- sample[s]
  sel.tree <- paste0(gsub("HSC_|hESC_","",sel.sample),".tree") %>% get()
  sel.bootstrap <- data.frame(sample=sel.sample, 
                              Bootstrap=sapply(str_split(sel.tree$node.label, "/"), "[", 2)) %>% na.omit()
  cat(sel.sample, "median", median(as.numeric(sel.bootstrap$Bootstrap)),"\n")
  #plot
  p <- ggplot(sel.bootstrap, aes(x = as.numeric(Bootstrap))) + geom_histogram(binwidth=1,fill="#1D91C0",alpha=0.5)+ #"#619EFF"
    labs(title=sel.sample, x="Bootstrap (%)",y="Counts")+theme_classic()+
    scale_x_continuous(breaks = seq(0,100,20))+
    theme(plot.title=element_text(size=12,color="black",hjust=0.5),axis.title=element_text(size=12,color="black"),
          axis.text=element_text(size=10,color='black'),
          axis.ticks=element_line(color="black"),axis.ticks.length=unit(0.13,"cm"))+
    geom_vline(xintercept=median(as.numeric(sel.bootstrap$Bootstrap)), color="#D73027", linetype = "dashed")
  #geom_text(x = 80, y = 20, label = "median", color = "black",family = "Times New Roman")
  return(p)
})

(p.bootstrap[[4]]/p.bootstrap[[5]]/p.bootstrap[[6]]/p.bootstrap[[7]])|
  (p.bootstrap[[1]]/p.bootstrap[[2]]/p.bootstrap[[3]]/p.bootstrap[[3]])

