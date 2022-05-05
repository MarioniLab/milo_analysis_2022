suppressPackageStartupMessages({
    library(miloR)
    library(SingleCellExperiment)
    library(Matrix)
    library(scater)
    library(readr)
    library(dplyr)
    library(patchwork)
    library(dyntoy)
    library(tidyr)
    library(tibble)
    library(BiocNeighbors)
    library(matrixStats)
    library(igraph)
    library(profmem)
})

print("importing simulation")
simulation <- readRDS("/nfs/research/marioni/mdmorgan/milo_testing/data/Merge_simulation_Milo.RDS")

## Make milo object
sim_milo <- Milo(simulation)
#print(names(colData(sim_milo)))

reducedDim(sim_milo, "PCA") <- reducedDim(sim_milo, "PCA")

reduced_dim <- rep(NA, 4)
graph <- rep(NA, 4)
n_cells <- rep(NA, 4)
j <- 1

for (i in c(190000, 142500, 95000, 47500)) {

     #down sample the milo object
     traj_milo <- sim_milo[,sample(sim_milo@colData@rownames, i)]

     #create design matrix to explain experimental design
     traj_design <- data.frame(colData(traj_milo))[,c("Sample", "Condition")]
     traj_design <- distinct(traj_design)
     rownames(traj_design) <- traj_design$Sample
     
     ### using system.time() calculate elapsed time
     a <- system.time({traj_milo4 <-buildGraph(traj_milo, k = 10, d = 30);
     traj_milo4 <- makeNhoods(traj_milo4, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme = "reduced_dim");
     traj_milo4 <- countCells(traj_milo4, meta.data = data.frame(colData(traj_milo)), sample="Sample");
     traj_milo4 <- calcNhoodDistance(traj_milo4, d=30);
     traj_milo4 <- testNhoods(traj_milo4, design = ~ Condition, design.df = traj_design, fdr.weighting = "k-distance")})

     b <- system.time({traj_milo5 <- buildGraph(traj_milo, k = 10, d = 30);
     traj_milo5 <- makeNhoods(traj_milo5, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme = "graph");
     traj_milo5 <- countCells(traj_milo5, meta.data = data.frame(colData(traj_milo)), sample="Sample");
     traj_milo5 <- testNhoods(traj_milo5, design = ~ Condition, design.df = traj_design, fdr.weighting = "graph-overlap")})
     
     n_cells[j] <- i
     reduced_dim[j] <- a[3]
     graph[j] <- b[3]
     j <- j + 1
 }
 
time.df <- cbind(n_cells, reduced_dim, graph)
time.df <- data.frame(time.df)
write_csv(time.df, "/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_2/time_df.csv")

# same as above, but for different values of k
# reduced_dimk <- rep(NA, 4)
# graphk <- rep(NA, 4)
# k_val <- rep(NA, 4)
# j <- 1
# 
# # create design matrix to explain experimental design
# traj_milo <- sim_milo[,sample(sim_milo@colData@rownames, 100000)]
# traj_design <- data.frame(colData(traj_milo))[,c("Sample", "Condition")]
# traj_design <- distinct(traj_design)
# rownames(traj_design) <- traj_design$Sample
# print(traj_design)
# 
# for (k in seq(10, 70, 20)) {
# 
#     print(k)
#     
#     ### using system.time() calculate elapsed time
#     c <- system.time({traj_milo2 <- miloR::buildGraph(traj_milo, k = k, d = 30);
#     traj_milo2 <- makeNhoods2(traj_milo2, prop = 0.1, k = k, d=30, refined = TRUE, refinement_scheme = "reduced_dim");
#     traj_milo2 <- miloR::countCells(traj_milo2, meta.data = data.frame(colData(traj_milo2)), sample="Sample");
#     traj_milo2 <- miloR::calcNhoodDistance(traj_milo2, d=30);
#     traj_milo2 <- testNhoods2(traj_milo2, design = ~ Condition, design.df = traj_design, fdr.weighting = "k-distance")})
# 
#     d <- system.time({traj_milo1 <- miloR::buildGraph(traj_milo, k = k, d = 30);
#     traj_milo1 <- makeNhoods2(traj_milo1, prop = 0.1, k = k, d=30, refined = TRUE, refinement_scheme = "graph");
#     traj_milo1 <- miloR::countCells(traj_milo1, meta.data = data.frame(colData(traj_milo1)), sample="Sample");
#     traj_milo1 <- testNhoods2(traj_milo1, design = ~ Condition, design.df = traj_design, fdr.weighting = "graph-overlap")})
# 
#     k_val[j] <- k
#     reduced_dimk[j] <- c[3]
#     graphk[j] <- d[3]
#     j <- j + 1
# }
# 
# k.df <- cbind(k_val, reduced_dimk, graphk)
# k.df <- data.frame(k.df)
# write_csv(k.df, "/nfs/research/marioni/akluzer/k_df.csv")