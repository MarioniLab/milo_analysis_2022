library(miloR)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix)
library(scater)
library(dplyr)
library(patchwork)
library(tidyr)
library(tibble)
library(igraph)
library(readr)
library(matrixStats)
library(BiocNeighbors)
library(scran)
library(miloR)

million_cells <- readRDS("/nfs/research/marioni/alsu/hubmap_metaRef/data/lung/sce_mapped.Rds")
print("million cells imported")

million_cells <- million_cells[,apply(reducedDim(million_cells, "pca.corrected"), 1, function(x) !all(is.na(x)))]

#reformat to create testable milo object
million_cells <- Milo(million_cells)
million_cells@colData$sample2 <- paste(million_cells@colData$dataset_origin, million_cells@colData$health_status, sep = "_")
unique(million_cells@colData$sample2)

## make design matrix
design <- data.frame(colData(million_cells))[,c("sample2", "health_status")]
design <- distinct(design)
rownames(design) <- design$sample2

# time for different values of k
graph <- rep(NA, 3)
k_val <- rep(NA, 3)
j <- 1

for (k in seq(10, 50, 20)) { 
    print(k)
    ### using system.time() calculate elapsed time
    a <- system.time({million_cells2 <- buildGraph(million_cells, k = k, d = 30, reduced.dim = "pca.corrected");
                                million_cells3 <- makeNhoods(million_cells2, prop = 0.1, k = k, d=30, refined = TRUE, reduced_dims = "pca.corrected", refinement_scheme = "graph");
                                million_cells4 <- countCells(million_cells3, meta.data = as.data.frame(colData(million_cells3)), samples="sample2");
                                da_results <- testNhoods(million_cells4, reduced.dim= "pca.corrected", design = ~ health_status, design.df = design, fdr.weighting = "graph-overlap")})
    
    k_val[j] <- k
    graph[j] <- a[3]
    j <- j + 1
}

time.df <- cbind(k_val, graph)
time.df <- data.frame(time.df)
write_csv(time.df, "/nfs/research/marioni/akluzer/speed_benchmarks/time_millioncells.csv")