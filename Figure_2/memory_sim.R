suppressPackageStartupMessages({
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
    library(devtools)
    library(miloR)
})

print("importing simulation")
simulation <- readRDS("/nfs/research/marioni/mdmorgan/milo_testing/data/Merge_simulation_Milo.RDS")

### ANAYSIS ###

## Make milo object
sim_milo <- Milo(simulation)
#sim_milo
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
    
    #Rprof(filename = paste0("Rprof_red", i, ".out"), interval = 0.02, memory.profiling = TRUE)
    gc1 <- sum(gc(reset = TRUE, full = TRUE)[,6])
    traj_milo4 <- buildGraph(traj_milo, k = 10, d = 30)
    traj_milo4 <- makeNhoods(traj_milo4, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme = "reduced_dim")
    traj_milo4 <- countCells(traj_milo4, meta.data = data.frame(colData(traj_milo)), sample="Sample")
    traj_milo4 <- calcNhoodDistance(traj_milo4, d=30)
    traj_milo4 <- testNhoods(traj_milo4, design = ~ Condition, design.df = traj_design, fdr.weighting = "k-distance")
    gc2 <- sum(gc(reset = FALSE, full = TRUE)[,6])
    red <- gc2 - gc1
    #Rprof(NULL)
    
    #Rprof(filename = paste0("Rprof_graph", i, ".out"), interval = 0.02, memory.profiling = TRUE)
    gc3 <- sum(gc(reset = TRUE, full = TRUE)[,6])
    traj_milo5 <- buildGraph(traj_milo, k = 10, d = 30)
    traj_milo5 <- makeNhoods(traj_milo5, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme = "graph")
    traj_milo5 <- countCells(traj_milo5, meta.data = data.frame(colData(traj_milo)), sample="Sample")
    traj_milo5 <- testNhoods(traj_milo5, design = ~ Condition, design.df = traj_design, fdr.weighting = "graph-overlap")
    gc4 <- sum(gc(reset = FALSE, full = TRUE)[,6])
    graph <- gc4 - gc3
    #Rprof(NULL)
    
    n_cells[j] <- i
    reduced_dim[j] <- red
    graph[j] <- graph
    j <- j + 1
}

mem.df <- cbind(n_cells, reduced_dim, graph)
mem.df <- data.frame(mem.df)
write_csv(mem.df, "/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_2/mem.df")