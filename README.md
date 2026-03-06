# PASTRI - Phylogenetic Adjacency-based State Transition Rate Inference

PASTRI is a computational framework designed to quantify cell state transition rates from phylogenetic data.  PASTRI utilizes cell phylogeny where terminal nodes (tips) are annotated by single-cell phenotypes, thereby inferring a transition rate matrix to robustly and systematically characterize cell-state transition dynamics across diverse biological contexts.

## System requirement
* Dependent packages: ape, argparse, data.table, expm, ggplot2, ggridges, ggtree, magrittr, Matrix, parallel, patchwork, pracma, reshape2, stats, tidyr,  tidyverse
* Require R (>= 4.2.3).
   
## Install
```
install.packages("devtools")
devtools::install_github("yuxiaochen11/PASTRI")
```

## Quickstart
```
library("PASTRI")
```
- The following files are needed for PASTRI.
1. A tree file of class "phylo" with node labels.
> ((((Cell_1:1,(Cell_2:1,Cell_3:1)1:1)1:1,Cell_4:1)1:1,(Cell_5:1,Cell_6:1,Cell_7:1)1:1)1:1,(Cell_8:1,Cell_9:1)1:1);
2. A dataframe with columns *TipLabel* and *TipAnn*, representing tip labels on the tree file and corresponding cell annotations.

> | nodeLabel | cell_ID | celltype | cellNum |
> | --- | --- | --- | --- |
> | Cell_1 | Cell_1 | C1 | 1 |
> | Cell_2 | Cell_2 | C1 | 2 |
> | Cell_2 | Cell_10 | C2 | 2 |
> | Cell_3 | Cell_3 | C2 | 1 |
> | Cell_4 | Cell_4 | C3 | 1 |
> | Cell_5 | Cell_5 | C0 | 1 |
> | Cell_6 | Cell_6 | C1 | 2 |
> | Cell_6 | Cell_11 | C1 | 2 |
> | Cell_7 | Cell_7 | C0 | 1 |
> | Cell_8 | Cell_8 | C1 | 1 |
> | Cell_9 | Cell_9 | C0 | 1 |

### Estimate transition rate with exemplar dataset.
- Get the depth fetures of tip (cell) pairs.
Using the function calculate_lca_depths() to compute pairwise LCA depth / height features from the specified tree file (sample.nwk). 
```
tree_file <- system.file("extdata", "HSC.nwk", package = "PASTRI")
node_depth_df <- calculate_lca_depths(tree.path = tree_file)
```
- Infer optimal constrained transition matrix.
Using the function get_optimal_transition_matrix() to obtain the optimal transition matrix. We will illustrate multiple iterations over a set of depth thresholds (depth) and combine all results into a single data frame. 
> **** 1. Read the CSV file that contains cell information.
```
cell_file <- system.file("extdata", "HSC_nodeInfos.csv", package = "PASTRI")
cell_info  <- read.csv(cell_file, stringsAsFactors = FALSE)
```
> **** 2. Specify a series of depth thresholds to test.
```
depth <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
```
> **** 3. Set cell type.
```
cell_type <- c("C4", "C5", "C6", "C7", "C8")
```
> **** 4. Constrained matrix. Represent prior knowledge such as irreversible state transitions.
```
Bound_matrix <- matrix(Inf, nrow = length(cell_type), ncol = length(cell_type))
Bound_matrix[upper.tri(Bound_matrix)] <- 0
```
> **** 5. Infer optimal constrained transition matrix.
```
df_depth <- lapply(seq(length(depth)), function(d) {
  sel.depth <- depth[d]
  print(sel.depth)  # Print the current depth threshold
    
  # Call get_optimal_transition_matrix() with the specified parameters
  res <- get_optimal_transition_matrix(
    node_pair_depth = node_depth_df,
    cell_info = cell_info,
    Sel_u = "lca_normalized_height",
    fi_depth = sel.depth,
    Bound_Matrix = Bound_matrix,
    cell_type_list = cell_type,
    mc.cores = 80
    )
  return(res)
})
```
> **** 6. Show the transition matrix of depth 0.5.
```
df_depth[[4]]$optimal_matrix %>% t()
```
## Contributing
### Contributors
Xiaochen Yu, https://yuxiaochen11.github.io/PASTRI/.

## Citations
When using PASTRI please cite:
- Yang W, Li Z, Yu X, et al. Inference of Stage-specific Transition Dynamics of Cellular States by Division History. 

## Note
PASTRI is still under development, some API names and call methods in notebooks may have changed, and the model performance may also be slightly different.
