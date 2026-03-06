rm(list=ls())
## ==================== import some librarys ==================================================
suppressMessages({
  library(patchwork)
  library(reshape2)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggbreak)
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(ggrastr)
  library(parallel)
  library(argparse)
  library(dplyr)
  library(magrittr)
  library(ggridges)
  library(stats)
  library(Matrix)
  library(expm)
  library(pracma)
  library(data.table)
  library(Seurat)
})

## ==================== source and define functions ===========================================
subtree.depths <- function(phylo.object){
  tree.tibble <- as_tibble(phylo.object) %>% as.data.frame()
  all.subtrees <- subtrees(phylo.object)
  distMatrix <- dist.nodes(phylo.object)
  # build in function to calculate depth
  getDist <- function(leaf.label, node){
    n2 <- dplyr::filter(tree.tibble, label == leaf.label) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  # get all subtrees depth
  all.subtrees.depth.df <- lapply(seq(length(all.subtrees)), function(s){
    subtree <- all.subtrees[[s]]
    subtree.root <- subtree$name
    subtree.depth <- Map(getDist, subtree$tip.label, subtree$name) %>% unlist() %>% as.vector() %>% max()
    return(data.frame(subtree.root = subtree.root,
                      depth = subtree.depth))
  }) %>% bind_rows()
  # return results
  return(all.subtrees.depth.df)
}

get.adistVStdist <- function(allele.info_path, tree_path){
  allele.info <- read.csv(allele.info_path, stringsAsFactors = F) %>% select(nodeLabel, editMutation) %>% unique()
  rownames(allele.info) <-  allele.info$nodeLabel
  node.edits <- allele.info %>% dplyr::select(-nodeLabel) %>% as.matrix()
  
  # get node compare df
  node.comp.df <- t(combn(rownames(node.edits), 2)) %>% as.data.frame(., stringsAsFactors =F) 
  colnames(node.comp.df) <- c("node1", "node2")
  # get tree depth as distance
  tree.obj <- read.tree(tree_path)
  distMatrix <- dist.nodes(tree.obj)
  treeTibble <- as_tibble(tree.obj) %>% as.data.frame()
  #all.subtree.depth <- subtree.depths(tree.obj)
  dist.comp <- 
    mclapply(seq(nrow(node.comp.df)), 
             mc.cores = 60,
             function(n){
               #n <- 1
               n1 <- node.comp.df$node1[n]
               n2 <- node.comp.df$node2[n]
               # extract node1 and node2 edits
               n1.edit <- strsplit(node.edits[n1,] %>% as.vector(), ",")[[1]]
               n2.edit <- strsplit(node.edits[n2,] %>% as.vector(), ",")[[1]]
               # compare allelic distance between these two nodes
               diff.num <- length(setdiff(n1.edit, n2.edit)) + length(setdiff(n2.edit, n1.edit))
               # compare tree dist (use the MRCA's depth as distance for two node)
               mrca <- getMRCA(tree.obj, c(n1, n2))
               #mrca.depth <- filter(all.subtree.depth, subtree.root == mrca) %>% `[`("depth") %>% unlist()
               #tree dist used the depth of max distrance of a cell pair to the most recent common ancestors
               #(the number of internal nodes on the path from one cell to the other)
               n1.seqName <- filter(treeTibble, label == n1) %>% `[`("node") %>% unlist()
               n2.seqName <- filter(treeTibble, label == n2) %>% `[`("node") %>% unlist()
               treedist <- max(distMatrix[n1.seqName, mrca], distMatrix[n2.seqName, mrca])
               # return results
               return(data.frame(node1 = n1,
                                 node2 = n2,
                                 allelic.dist = diff.num,
                                 tree.dist = treedist,
                                 stringsAsFactors = F))
             }) %>% bind_rows()
  
  dist.comp$nor.allelic.dist <- dist.comp$allelic.dist / max(dist.comp$allelic.dist)
  dist.comp$nor.tree.dist <- dist.comp$tree.dist / max(dist.comp$tree.dist)
  dist.comp$node1 <- as.character(dist.comp$node1)
  dist.comp$node2 <- as.character(dist.comp$node2)
  return(dist.comp)
}

checkCon2 <- function(Node1, Node2, alleleInfo, rawEvents){
  con1 <- filter(alleleInfo, nodeLabel == Node1)
  con2 <- filter(alleleInfo, nodeLabel == Node2)
  # get cell BC
  BC1 <- con1$nodeLabel[1]
  BC2 <- con2$nodeLabel[1]
  # get raw events
  raw1 <- filter(rawEvents, BC == BC1) %>% `[`("editMutation") 
  raw2 <- filter(rawEvents, BC == BC2) %>% `[`("editMutation") 
  # overlap edits
  overlap <- dplyr::intersect(raw1, raw2)
  return(data.frame(node1 = Node1,
                    node2 = Node2,
                    overlapE = nrow(overlap), stringsAsFactors = F))
}

getOverlapInfo.dist <- function(Pair.df, alleleInfo, rawEvents){
  OverlapInfo <- mclapply(seq(nrow(Pair.df)), function(i){
    # i <- 1
    n1 <- Pair.df$node1[i]
    n2 <- Pair.df$node2[i]
    treeD <- Pair.df$treeDist[i]
    subCheck <- checkCon2(Node1 = n1,
                          Node2 = n2,
                          alleleInfo = alleleInfo,
                          rawEvents = rawEvents)
    subCheck$treeDist <- treeD
    return(subCheck)
  }, mc.cores = 60) %>% bind_rows()
  return(OverlapInfo)
}

run_get_results <- function(alleleInfos_path, moreInfos_path, oneCellNode.pair.dist_path, checkString){
  # input file
  alleleInfos <- read.csv(alleleInfos_path, stringsAsFactors = F)
  moreInfos <- read.csv( moreInfos_path, sep = '\t', stringsAsFactors = F)
  moreInfos$BC <- gsub("=", "_", moreInfos$BC)
  #oneCellNode.pair.dist <- read.csv(oneCellNode.pair.dist_path, stringsAsFactors = F)
  oneCellNode.pair.dist <- readRDS(oneCellNode.pair.dist_path) %>% select(node1, node2, nor.tree.dist)
  colnames(oneCellNode.pair.dist) <- c("node1","node2","treeDist")
  
  # check have anc cells and oneCell node
  haveANC.cells <- moreInfos %>% filter(., (!!sym(checkString)) == "anc") %>% `[`("BC") %>% unlist() %>% unique() # 
  haveANC.oneCellNode <- alleleInfos %>% filter(nodeLabel %in% haveANC.cells & cellNum == 1) %>% `[`("nodeLabel") %>% unlist() # 
  
  # check overlap infos
  moreInfos.no.dis <- moreInfos %>% filter(., (!!sym(checkString)) != "dis")
  overlapInfo.dist <- getOverlapInfo.dist(Pair.df = oneCellNode.pair.dist,
                                          alleleInfo = alleleInfos,
                                          rawEvents = moreInfos.no.dis)
  # check share_ratio
  if (length(haveANC.oneCellNode) == 1){
    haveANC.share_ratio <- NULL
  }else{
    #haveANC.pair.dist <- oneCellNode.pair.dist %>% filter(node1%in%haveANC.oneCellNode & node2%in%haveANC.oneCellNode)
    #haveANC.overlapInfo <- getOverlapInfo.dist(Pair.df = haveANC.pair.dist,
    #                                           alleleInfo = alleleInfos,
    #                                           rawEvents = moreInfos.no.dis)
    haveANC.overlapInfo <- overlapInfo.dist %>% filter(node1%in%haveANC.oneCellNode & node2%in%haveANC.oneCellNode) %>%
      mutate(treeDist = round(treeDist, 8))
    haveANC.share_ratio <- haveANC.overlapInfo %>% dplyr::group_by(treeDist) %>% 
      dplyr::summarise(share_ratio=sum(overlapE > 0) / n()) %>% ungroup()
  }
  
  return(list(haveANC_cells=haveANC.cells,
              haveANC_oneCellNode=haveANC.oneCellNode,
              Cells=nrow(alleleInfos),
              Nodes=length(unique(alleleInfos$nodeLabel)),
              oneCellNodes=nrow(filter(alleleInfos, cellNum == 1)),
              overlapInfo_dist=overlapInfo.dist,
              haveANC_share_ratio=haveANC.share_ratio))
}

get.latest.common.ancestor.depth <- function(tree.path,allele.info_path){
  path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/"
  tree <- read.tree(paste0(path, tree.path))
  #==subtree.node to root height
  df_tree <- fortify(tree)
  all_parents <- unique(df_tree$parent)
  subtree.depth <- mclapply(1:length(all_parents),function(j){
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
  })%>% plyr::rbind.fill()
  
  subtree.depth$Dt <- max(subtree.depth$Ds)
  subtree.depth$Ds.norm <- subtree.depth$Ds/subtree.depth$Dt
  subtree.depth$mdepth <- ((subtree.depth$Dr/subtree.depth$Dt)+(1-subtree.depth$Ds/subtree.depth$Dt))/2
  subtree.depth$node.height <- 1-subtree.depth$mdepth
  
  #==calcultate single leaf (endpoint cell) pairs depth
  print("==== calculate node pair depth and normalized height ====")
  allele.info <- read.csv(paste0(path, allele.info_path), stringsAsFactors = F) %>% select(nodeLabel,BC)
  node.node_pair_com_df <- combn(unique(allele.info$nodeLabel %>% unlist()), 2)
  node.comp.df <- data.frame(node1 = node.node_pair_com_df[1,],
                             node2 = node.node_pair_com_df[2,], 
                             stringsAsFactors =F)
  #node.comp.df <- t(combn(unique(allele.info$nodeLabel), 2)) %>% as.data.frame(., stringsAsFactors =F)
  
  dist.comp <- mclapply(seq(nrow(node.comp.df)), mc.cores = 60, function(n){
    n1 <- node.comp.df$node1[n]
    n2 <- node.comp.df$node2[n]
    #latest common ancestor depth (use the MRCA's depth as distance for two node)
    mrca <- getMRCA(tree, c(n1, n2))
    mrca.depth <- filter(subtree.depth, subtree.root == mrca) %>% `[`("Ds") %>% unlist()
    mrca.depth.norm <- filter(subtree.depth, subtree.root == mrca) %>% `[`("Ds.norm") %>% unlist()
    mrca.mheight <- filter(subtree.depth, subtree.root == mrca) %>% `[`("node.height") %>% unlist()
    # return results
    return(data.frame(node1 = n1,
                      node2 = n2,
                      subtree.root = mrca,
                      mrca.depth = mrca.depth,
                      mrca.depth.norm = mrca.depth.norm,
                      mrca.mheight = mrca.mheight,
                      stringsAsFactors = F))
  }) %>% bind_rows()
  dist.comp$node1 <- as.character(dist.comp$node1)
  dist.comp$node2 <- as.character(dist.comp$node2)
  return(dist.comp)
}

get.tree.expMean.dist <- function(treeDist,allele.info_path,exp_path){
  #get onecell node pair treeDist
  path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/"
  allele.info <- read.csv(paste0(path, allele.info_path), stringsAsFactors = F)
  onecell.node <- filter(allele.info,allele.info$cellNum==1)
  onecell.treeDist <- treeDist[treeDist$node1%in%onecell.node$nodeLabel&treeDist$node2%in%onecell.node$nodeLabel,]
  #get exp matrix, filter more than 50% cell exp genes or hvg.genes
  exp.matrix <- fread(exp_path,sep=',',header=T) %>% as.data.frame()
  rownames(exp.matrix) <- exp.matrix$BC
  exp.matrix <- exp.matrix %>% select(-c(V1,BC)) %>% t() %>% as.data.frame()
  exp.matrix$freq <- apply(exp.matrix,1,function(x){sum(x>0)/ncol(exp.matrix)})
  exp.matrix <- exp.matrix[exp.matrix$freq>0.5,] %>% select(-ncol(exp.matrix))
  #exp.matrix <- exp.matrix %>% select(-ncol(exp.matrix))
  #exp.matrix <- exp.matrix[rownames(exp.matrix) %in% hvg.genes$V1,]
  
  treeDistAndExpDist <- mclapply(seq(nrow(onecell.treeDist)), function(j){
    node1 <- onecell.treeDist$node1[j]
    node2 <- onecell.treeDist$node2[j]
    treeDist1 <- onecell.treeDist$mrca.depth[j]
    treeDist2 <- onecell.treeDist$mrca.mheight[j]
    subtree.root <- onecell.treeDist$subtree.root[j]
    # get exp 
    node1.exp <- exp.matrix[, node1] %>% as.numeric()
    node2.exp <- exp.matrix[, node2] %>% as.numeric()
    #expDist <- dist(rbind(node1.exp, node2.exp), method = "euclidean")[1]
    #expDist2 <- 1-cor(node1.exp, node2.exp,method="pearson")
    expDist <- amap::Dist(rbind(node1.exp,node2.exp),method="euclidean")[1]
    expDist2 <- amap::Dist(rbind(node1.exp,node2.exp),method="pearson")[1]
    return(data.frame(node1=node1,
                      node2=node2,
                      subtree.root=subtree.root,
                      treeDist=treeDist1,
                      mdepth=treeDist2,
                      expDist.eu=expDist,
                      expDist.cor=expDist2,
                      stringsAsFactors = F))
  }, mc.cores = 60) %>% bind_rows(.)
  
  #sampled_eu <- mean(treeDistAndExpDist$expDist.eu[sample(nrow(treeDistAndExpDist), 1000, replace = FALSE)])
  #sampled_cor <- mean(treeDistAndExpDist$expDist.cor[sample(nrow(treeDistAndExpDist), 1000, replace = FALSE)])
  sub.node <- unique(treeDistAndExpDist$subtree.root)
  final.E.L.dist <- mclapply(seq(length(sub.node)),function(x){
    select.sub.node <- sub.node[x]
    local.subtree <- treeDistAndExpDist[treeDistAndExpDist$subtree.root==select.sub.node,]
    local.subtree.dep <- as.numeric(unique(local.subtree$treeDist))
    local.subtree.mdep <- as.numeric(unique(local.subtree$mdepth))
    local.ExpDist.mean <- data.frame(subtree.root=select.sub.node,
                                     depth = local.subtree.dep,
                                     mdepth = local.subtree.mdep,
                                     expDist.eu=mean(local.subtree$expDist.eu),
                                     expDist.cor=mean(local.subtree$expDist.cor),
                                     stringsAsFactors = F) #%>% mutate(norm.eu=expDist.eu/sampled_eu,norm.cor=expDist.cor/sampled_cor)
  }) %>% bind_rows(.)
  return(final.E.L.dist)
}


## ============Run=============================================================================
#===================================================================================
#hESC_1 (hESC_A7), hESC_2 (hESC_C9), hESC_3 (hESC_C12), hESC_4 (hESC_D3)
#HSC_1 (HSC_B7), HSC_2 (HSC_E12), HSC_3 (HSC_F12)
#===================================================================================

##==Figure 4H and Figure S6E, compare allelicDist with treeDist
input_path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/"
#get alle dist and tree dist
A7.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "hESC_A7_B/hESC_A7.AllelesInfo.csv"),
                                    tree_path=paste0(input_path, "hESC_A7_B/hESC_A7.nwk")) %>% mutate(sample="hESC_A7")
C9.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "hESC_C9_B/hESC_C9.AllelesInfo.csv"),
                                    tree_path=paste0(input_path, "hESC_C9_B/hESC_C9.nwk")) %>% mutate(sample="hESC_C9")
C12.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "hESC_C12_B/hESC_C12.AllelesInfo.csv"),
                                     tree_path=paste0(input_path, "hESC_C12_B/hESC_C12.nwk")) %>% mutate(sample="hESC_C12")
D3.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "hESC_D3_B/hESC_D3.AllelesInfo.csv"),
                                    tree_path=paste0(input_path, "hESC_D3_B/hESC_D3.nwk")) %>% mutate(sample="hESC_D3")
B7.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "HSC_B7_B/HSC_B7.AllelesInfo.csv"),
                                    tree_path=paste0(input_path, "HSC_B7_B/HSC_B7.nwk")) %>% mutate(sample="HSC_B7")
E12.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "HSC_E12_B/HSC_E12.AllelesInfo.csv"),
                                     tree_path=paste0(input_path, "HSC_E12_B/HSC_E12.nwk")) %>% mutate(sample="HSC_E12")
F12.adistVStdist <- get.adistVStdist(allele.info_path=paste0(input_path, "HSC_F12_B/HSC_F12.AllelesInfo.csv"),
                                     tree_path=paste0(input_path, "HSC_F12_B/HSC_F12.nwk")) %>% mutate(sample="HSC_F12")
save(A7.adistVStdist,C9.adistVStdist,C12.adistVStdist,D3.adistVStdist,B7.adistVStdist,E12.adistVStdist,F12.adistVStdist,
     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/4.0.1.AllelicDist_TreeDist.Rda")

#compare alle dist and tree dist
rm(list=ls())
load("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/4.0.1.AllelicDist_TreeDist.Rda")
sample <- c("hESC_A7","hESC_C9","hESC_C12","hESC_D3","HSC_B7","HSC_E12","HSC_F12")

res.adistVStdist <- lapply(seq(sample), function(s){
  #s <- 1
  sel.sample <- sample[s]
  sel.adistVStdist <- paste0(sub("^.*_", "", sel.sample), ".adistVStdist") %>% get()
  res.rho <- cor.test(sel.adistVStdist$nor.allelic.dist,sel.adistVStdist$nor.tree.dist,method = "spearman")
  res.cor <- cor.test(sel.adistVStdist$nor.allelic.dist,sel.adistVStdist$nor.tree.dist,method = "pearson")
  
  res <- data.frame(sample=sel.sample, stringsAsFactors =F,
                    rho=res.rho$estimate, rho.pval=res.rho$p.value,
                    cor=res.cor$estimate, cor.pval=res.cor$p.value)
  return(res)
}) %>% bind_rows()
res.adistVStdist$sample <- factor(res.adistVStdist$sample, levels = sample)

#res.adistVStdist.new <- res.adistVStdist %>% select(sample, rho, cor) %>% melt(id="sample")
#res.adistVStdist.new$sample <- factor(res.adistVStdist.new$sample, levels=sample)
p.rho <- ggplot(res.adistVStdist, aes(x=sample, y=rho)) + 
  geom_bar(stat="identity",fill="#1D91C0",alpha=0.5,color="black",width=0.6)+
  labs(x="",y="Correlation between lineage distance \nand allelic distance of barcode")+theme_classic()+scale_y_continuous(expand=c(0,0))+
  theme(axis.title=element_text(size=12,color="black"),axis.text.y=element_text(size=12,color='black'),
        axis.text.x=element_text(size=12,color='black', angle=30, hjust=1))
p.rho

p.bin2d <- lapply(seq(sample), function(s){
  #s <- 1
  sel.sample <- sample[s]
  sel.adistVStdist <- paste0(sub("^.*_", "", sel.sample), ".adistVStdist") %>% get()
  res.rho <- cor.test(sel.adistVStdist$nor.allelic.dist,sel.adistVStdist$nor.tree.dist,method = "spearman")

  #
  p <- ggplot(sel.adistVStdist, aes(x=nor.tree.dist, y=nor.allelic.dist))+
    geom_bin2d(aes(fill=after_stat(density)), bins = 10)+scale_fill_gradient(low="white",high="red")+
    labs(title = paste0("\nSpearman rho = ", round(res.rho$estimate,2), " P < 10e-300"),
         x="Normalized lineage distance", y="Normalized allelic distance")+
    scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
    theme(panel.border=element_rect(fill='transparent', color='black'), panel.background=element_blank(),
          plot.title=element_text(size=10,color="black",), axis.title=element_text(size=10,color="black"),
          axis.text=element_text(size=10,color="black"), axis.line=element_line(color="black"),
          legend.title=element_text(size=10), legend.text=element_text(size=10))
  return(p)
})
(p.bin2d[[1]]|p.bin2d[[2]]|p.bin2d[[3]]|p.bin2d[[4]])/(p.bin2d[[5]]|p.bin2d[[6]]|p.bin2d[[7]])

##==Figure 4I, ANC_share_ratio_uniqueEvents
#get ANC share ratio
sample <- c("hESC_A7","hESC_C9","hESC_C12","hESC_D3","HSC_B7","HSC_E12","HSC_F12")

res.anc.share <- lapply(seq(sample), function(s){
  #s <- 1
  sel.sample <- sample[s]
  print(sel.sample)
  input_path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/"
  
  sel.ANC.share.results <- run_get_results(
    alleleInfos_path = paste0(input_path, "6.iqtree/",sel.sample,"_B/",sel.sample,".AllelesInfo.csv"),
    moreInfos_path = paste0(input_path, "5.ConEvents/", sel.sample, "_moreInfos.txt"),
    oneCellNode.pair.dist_path = paste0(input_path, "6.iqtree/",sel.sample,"_B/",sel.sample,".oneCellNode.pair.dist.Rds"),
    checkString = "checkInfo")
  
  sel.anc.share <- sel.ANC.share.results$haveANC_share_ratio %>% mutate(sample=sel.sample,
                                                                        Type=sub("_.*", "", sel.sample))
  print("Success")
  return(sel.anc.share)
}) %>% bind_rows()


res.anc.share.new <- res.anc.share %>% mutate(nor.treeDist.group=case_when(treeDist<=0.2~0.2,
                                                                           treeDist<=0.4&treeDist>0.2~0.4,
                                                                           treeDist<=0.6&treeDist>0.4~0.6,
                                                                           treeDist<=0.8&treeDist>0.6~0.8,
                                                                           treeDist<=1.0&treeDist>0.8~1.0))
##hESC and HSC
res.anc.share.new2 <- res.anc.share.new %>% group_by(Type, nor.treeDist.group) %>% 
  dplyr::summarize(mean=mean(share_ratio),sd=sd(share_ratio),sd.err=sd(share_ratio)/sqrt(n())) %>% ungroup() %>% 
  as.data.frame() %>% mutate(Type=factor(Type, levels=c("hESC", "HSC")))

p.anc.share <- ggplot(res.anc.share.new2, aes(x=nor.treeDist.group, y=mean)) +
  geom_bar(stat="identity",color="black",fill="#CCEBC5",width = 0.1)+theme_classic()+
  geom_errorbar(aes(ymin=mean-sd.err, ymax=mean+sd.err), width=0.05)+
  facet_wrap(~Type, nrow=1)+#scale_y_break(c(0.095,0.1),scales = 0.3,space=0.2)+
  labs(x="Normalized lineage distance",y="Share ancestral Barcode ratio")+
  theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"),
        axis.line=element_line(color="black"))
p.anc.share


##==Figure_S6B, mutation type frequency
A7 <- read.table(paste0(input_path, "hESC_A7_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="hESC_A7",freq=Num/sum(Num),type="hESC")
C9 <- read.table(paste0(input_path, "hESC_C9_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="hESC_C9",freq=Num/sum(Num),type="hESC")
C12 <- read.table(paste0(input_path, "hESC_C12_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="hESC_C12",freq=Num/sum(Num),type="hESC")
D3 <- read.table(paste0(input_path, "hESC_D3_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="hESC_D3",freq=Num/sum(Num),type="hESC")
B7 <- read.table(paste0(input_path, "HSC_B7_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="HSC_B7",freq=Num/sum(Num),type="HSC")
E12 <- read.table(paste0(input_path, "HSC_E12_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="HSC_E12",freq=Num/sum(Num),type="HSC")
F12 <- read.table(paste0(input_path, "HSC_F12_MutType.txt"),header=T,sep = '\t') %>% mutate(sample="HSC_F12",freq=Num/sum(Num),type="HSC")

all <- rbind(A7, C9, C12, D3, B7, E12, F12) %>% mutate(MutType=gsub("count_", "", MutType)) %>%
  #group_by(type, MutType) %>% mutate(type_num=sum(Num)) %>% ungroup() %>% select(-c(Num, sample, freq)) %>%
  #unique() %>% group_by(type) %>% mutate(Freq=type_num/sum(type_num)) %>% ungroup()
  group_by(MutType) %>% mutate(type_num=sum(Num)) %>% ungroup() %>% select(MutType,type_num) %>% unique() %>%
  mutate(Freq=type_num/sum(type_num))
all$MutType <- factor(all$MutType,levels = unique(all$MutType))

p.MutType.freq <- ggplot(all, aes(x=MutType,y=Freq))+geom_bar(stat="identity", width=0.6, color="black", fill="gray")+
  theme_classic()+labs(x="Mutation type",y="Mutation frequency")+ #facet_grid(.~class)+
  scale_y_break(c(0.05,0.5),scales = 0.3,space=0.15,expand = c(0,0))+ #geom_col_pattern(pattern="forwardslash")+
  theme(axis.title=element_text(size=12,color="black"), panel.background=element_blank(),
        axis.text.x=element_text(size=10,color='black',angle=30,hjust=1,vjust=1),
        axis.text.y=element_text(size=10,color='black'),
        axis.ticks=element_line(color="black"),axis.ticks.length=unit(0.13,"cm"))
p.MutType.freq

##==Figure_S6C, mutation number per allele
path <- "/mnt/data5/disk/yangwj/4.1.Pacbio_blastr/4-6.filter_numIndel_zws/Indels20_zws2_7new/"
hESC_A7 <- read.table(paste0(path, "5.ConEvents/hESC_A7_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="hESC_A7")
hESC_C9 <- read.table(paste0(path, "5.ConEvents/hESC_C9_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="hESC_C9")
hESC_C12 <- read.table(paste0(path, "5.ConEvents/hESC_C12_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="hESC_C12")
hESC_D3 <- read.table(paste0(path, "5.ConEvents/hESC_D3_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="hESC_D3")
HSC_B7 <- read.table(paste0(path, "5.ConEvents/HSC_B7_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="HSC_B7")
HSC_E12 <- read.table(paste0(path, "5.ConEvents/HSC_E12_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="HSC_E12")
HSC_F12 <- read.table(paste0(path, "5.ConEvents/HSC_F12_EditconEvents.txt"),header=T,sep='\t')[,c(1,3)] %>% mutate(sample="HSC_F12")

data <- rbind(hESC_A7, hESC_C9, hESC_C12, hESC_D3, HSC_B7, HSC_E12,HSC_F12)
data$sample <- factor(data$sample,levels = unique(data$sample))
pMut_num <- ggplot(data, aes(x=sample,y=num)) + geom_violin(trim=TRUE,color="gray",fill="gray")+
  geom_boxplot(width=0.2,position=position_dodge(0.1),outlier.color = NA)+
  labs(x="",y="Mutation number per allele")+theme_classic()+ 
  #scale_y_continuous(limits = c(0,30),expand=c(0,0))+
  theme(panel.background=element_blank(),axis.ticks=element_line(color="black"),axis.ticks.length=unit(0.13,"cm"),
        axis.title.y = element_text(size=12,color="black"),axis.text.x=element_text(size=10,color='black',angle=30,hjust=1,vjust=1),
        axis.text.y=element_text(size=10,color="black"))
rm(data, hESC_A7, hESC_C9, hESC_C12, hESC_D3, HSC_B7, HSC_F12, HSC_E12)
pMut_num


##==Figure_S6F, compare transcriptomic difference with lineage distance
#example
#Get expression matrix per HSC sample
sce <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.ref.scobj.SOT.V4.Rds")
sce.ref.exp <- as.data.frame(t(as.matrix(sce@assays$RNA@data))) #logNorm标准化后数据，scale.data是缩放后的data
sce.ref.exp[1:10,1:4]
#B7_HSC
cell.exp <- sce.ref.exp[grep("^B7_HSC_", rownames(sce.ref.exp)),]
cell.exp <- cell.exp %>% mutate(BC=rownames(cell.exp),BC=gsub("B7_HSC_","BC_",BC),BC=gsub("-1","",BC))
a <- cell.exp %>% select(27312)
#The last column is BC (lineage barcode), and 1:27312 is the expression matrix
write.csv(cell.exp,"/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/expMatrix/B7_HSC_allcell.exp.csv")

#Get expression matrix per hESC sample
sce.query <- readRDS("/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/SOT/numBalance.query.scobj.SOT.V4.Rds")
hESC_ref.exp <- sce@assays$RNA@data %>% .[, grep("hESC", colnames(.))] %>% as.data.frame() %>% mutate(gene=rownames(.))
hESC_ref.exp[1:2,1:2] 
hESC_query.exp <- sce.query@assays$RNA@data %>% .[, grep("hESC", colnames(.))] %>% as.data.frame() %>% mutate(gene=rownames(.))
hESC_query.exp[1:2,1:2]
hESC.exp <- full_join(hESC_ref.exp, hESC_query.exp, by="gene") 
rownames(hESC.exp) <- hESC.exp$gene
hESC.exp <- hESC.exp %>% select(-gene) %>% t() %>% as.data.frame()
hESC.exp[1:2,1:2]
#A7_hESC
cell.exp <- hESC.exp[grep("^A7_hESC", rownames(hESC.exp)),]
cell.exp <- cell.exp %>% mutate(BC=rownames(cell.exp),BC=gsub("A7_hESC_","BC_",BC),BC=gsub("-1","",BC))
a <- cell.exp %>% select(27312)
write.csv(cell.exp,"/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/expMatrix/A7_hESC_allcell.exp.csv")

#Get lineage distance and transcriptomic difference
B7.node.depth <- get.latest.common.ancestor.depth(
  tree.path="4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7_B/HSC_B7.nwk",
  allele.info_path="4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7_B/HSC_B7.AllelesInfo.csv")
B7.E.L.dist <- get.tree.expMean.dist(
  treeDist=B7.node.depth,
  allele.info_path="4-6.filter_numIndel_zws/Indels20_zws2_7new/6.iqtree/HSC_B7_B/HSC_B7.AllelesInfo.csv",
  exp_path="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/expMatrix/B7_HSC_allcell.exp.csv")
#compare
cor.test(B7.E.L.dist$mdepth,B7.E.L.dist$expDist.cor,method="spearman")#rho 0.1861114 p-value = 0.0014

B7.E.L.dist <- B7.E.L.dist %>% mutate(mdepth.norm=case_when(
  mdepth<=0.2~0.2,mdepth<=0.3&mdepth>0.2~0.3,
  mdepth<=0.4&mdepth>0.3~0.4,mdepth<=0.5&mdepth>0.4~0.5,mdepth<=0.6&mdepth>0.5~0.6,mdepth<=0.7&mdepth>0.6~0.7,
  mdepth<=0.8&mdepth>0.7~0.8,mdepth<=0.9&mdepth>0.8~0.9,mdepth<=1.0&mdepth>0.9~1))
table(B7.E.L.dist$mdepth.norm)

ggplot(B7.E.L.dist,aes(x=mdepth.norm,y=expDist.cor,group=mdepth.norm))+geom_point(color="gray",size=1)+
  geom_boxplot(width=0.05,position=position_dodge(0.01),outlier.color = NA)+ylim(0,0.5)+
  labs(title = "B7_HSC rho 0.19 P < 10e-2",x="Normalized depth",y="Transcriptomic difference")+ theme_classic()+
  theme(plot.title=element_text(size=12, hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),
        axis.line=element_line(color="black"),axis.text=element_text(size=12,color="black"),axis.title=element_text(size=12,color="black"))

save(A7.node.depth, C9.node.depth, C12.node.depth, D3.node.depth, B7.node.depth, E12.node.depth, F12.node.depth,
     file="/mnt/data5/disk/yangwj/Result_RData_Rds/Cell_State_Transition/Fig4_HSC/hESC_HSC.Node.depth.Rda")

