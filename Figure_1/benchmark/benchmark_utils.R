### BENCHMARKING FUNCTIONS ###

suppressPackageStartupMessages({
  #BiocManager::install("DESeq")
  #devtools::install_github("MarioniLab/miloR", ref="devel", force = TRUE)
  library(SingleCellExperiment)
  library(DAseq)
  library(BiocNeighbors)
  library(Matrix)
  library(miloR)
  library(tibble)
  library(dplyr)
  library(igraph)
  library(cydar)
  library(pdist)
  library(reshape2)
  library(limma)
  library(edgeR)
  library(BiocGenerics)
})

# # ## Set-up reticulate 4 MELD
# reticulate::use_condaenv("emma_env", required=TRUE)
# library(reticulate) ## development version of reticulate, or numba use breaks C stack

### SYNTHETIC LABELS ###

.find_centroid <- function(X_emb, cluster_membership){
  cl.ixs <- split(1:nrow(X_emb), cluster_membership)  
  centroid_emb <- sapply(cl.ixs, function(x) colMeans(X_emb[x, , drop=FALSE]))
  centroid_emb
}

.member_weight <- function(x, centroid_dist, m=2){
  # centroid_dist <- pdist(t(x), t(centroid_emb))@dist
  w_memb <- sapply(centroid_dist, function(x) 1/sum(x/centroid_dist)^(2/(m-1)))
}

.scale_to_range <- function(x, min=1, max=10){
  ((x - min(x))/(max(x)-min(x)))*(max-min) + min
}


.logit <- function(x, a=1){
  1/(1+exp(- a * x))
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels_pop <- function(sce, # SingleCellExperiment obj
                                     pop, pop_column="celltype",
                                     pop_enr = 0.7,
                                     redDim='pca.corrected', # embedding to use to simulate differential abundance
                                     n_conditions=2, # number of conditions to simulate
                                     n_replicates=3, # number of replicates per condition
                                     n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                     condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                     m=2, # Fuzziness parameter (higher m, more fuzziness)
                                     a_logit=0.5, # logit parameter
                                     cap_enr=NULL,
                                     seed=42){
  
  # pop_sce = sce[,sce[[pop_column]]==pop]
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  X_emb = reducedDim(sce, redDim)
  
  ## Find cluster center
  cluster_membership = sce[[pop_column]]
  centroid_emb <- .find_centroid(X_emb, cluster_membership)
  
  ## Assign weight to each cell for each cluster center
  centroid_dist <- pdist(X_emb, t(centroid_emb))
  centroid_dist <- as.matrix(centroid_dist)
  
  w <- sapply(1:ncol(centroid_dist),  function(j) 
    sapply(1:nrow(centroid_dist), function(i) 
      1/sum(centroid_dist[i,j]/centroid_dist[i,])^(2/(m-1))
    ) 
  )
  colnames(w) <- colnames(centroid_emb)
  rownames(w) <- rownames(X_emb)
  w <- apply(scale(w), 2, .logit, a=a_logit)
  ## Normalize weights from enr_score to 0.5
  enr_scores <- rep(0.5, ncol(w)) ## Generate enrichment prob for each cluster
  # enr_scores <- runif(ncol(w)) ## Generate enrichment prob for each cluster
  names(enr_scores) <- colnames(w)
  if(length(pop_enr) == length(pop)){
    enr_scores[pop] <- pop_enr
  } else{
    # assume all pops have the same enrichment
    pop_enr <- rep(pop_enr, length(pop))
    enr_scores[pop] <- pop_enr
  }
  
  
  # altering the baseline probability can induce a skew towards a condition across _all_ cells
  enr_prob <- sapply(1:ncol(w), function(i) .scale_to_range(w[,i], min=0.5*condition_balance,
                                                            max=enr_scores[i]))
  colnames(enr_prob) <- colnames(centroid_emb)
  
  # need to integrate over these to get the condition probabilities
  # need to set relevant pops only, force the others to ~0.5
  prob_matrix <- enr_prob[,pop]
  if(is(prob_matrix, "matrix")){
    cond_probability <- rowMeans(prob_matrix)
    for(x in seq_along(pop)){
      cond_probability[sce[[pop_column]] == pop[x]] <- prob_matrix[sce[[pop_column]] == pop[x], pop[x]]
    }
  } else{
    cond_probability <- prob_matrix
  }
  
  ## Cap probabilities (to have the same number of DA cells w different maximum Fold Change)
  if (!is.null(cap_enr)) {
    cond_probability <- ifelse(cond_probability > cap_enr, cap_enr, cond_probability)
  }
  # sim3$Condition1_prob <- ifelse(sim3$Condition1_prob > 0.8, 0.8, sim3$Condition1_prob)
  
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
    names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

## Simple synthetic condition labelling based on cluster membership
add_synthetic_labels_by_cluster <- function(sce, # SingleCellExperiment obj
                                            pop, pop_column="celltype",
                                            pop_enr = 0.7,
                                            # redDim='pca.corrected', # embedding to use to simulate differential abundance
                                            n_conditions=2, # number of conditions to simulate
                                            n_replicates=3, # number of replicates per condition
                                            n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                            condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                            m=2, # Fuzziness parameter (higher m, more fuzziness)
                                            a_logit=0.5, # logit parameter
                                            cap_enr=NULL,
                                            seed=42){
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  ## Set prop != 0.5 for cells in pop
  cond_probability <- rep(0.5, ncol(sce))
  cond_probability[sce[[pop_column]] == pop] <- pop_enr
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
    names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

## Group true DA cells in clusters (to test coverage of DA methods)
cluster_synthetic_labels <- function(embryo_sce, graph, min_cl_size=5){
  adj <- get.adjacency(graph)
  
  ## Retain DA cells
  da.adj <- adj[embryo_sce$true_labels!='NotDA',embryo_sce$true_labels!='NotDA']
  
  ## REmove edges between cells with discodant LFC sign
  da.cells.mat <- sapply(unique(embryo_sce$true_labels), function(x) as.numeric(embryo_sce$true_labels==x))
  da.cells.mat <- da.cells.mat[embryo_sce$true_labels!='NotDA',c("NegLFC", "PosLFC")]
  concord.sign.adj <- tcrossprod(da.cells.mat[,c("NegLFC", "PosLFC")], da.cells.mat[,c("NegLFC", "PosLFC")])
  concord.sign.adj <- as(concord.sign.adj, 'sparseMatrix')
  da.adj[concord.sign.adj==0] <- 0
  
  ## Cluster DA cells
  da.graph <- graph_from_adjacency_matrix(da.adj, mode = 'undirected')
  clust <- igraph::cluster_louvain(da.graph)
  embryo_sce$true_DA_clust <- rep(NA, length(embryo_sce$true_labels))
  embryo_sce$true_DA_clust[embryo_sce$true_labels != "NotDA"] <- clust$membership
  
  ## Remove singletons (or less than min_cl_size cells)
  embryo_sce$true_DA_clust[embryo_sce$true_DA_clust %in% which(table(clust$membership) < min_cl_size)] <- NA
  
  embryo_sce
}

### SYNTHETIC BATCH EFFECT ###

add_batch_effect <- function(embryo_sce, batch_col="synth_samples", norm_sd=0.5){
  cellids_sample <- split(embryo_sce$cell, embryo_sce[[batch_col]])
  X_pca <- reducedDim(embryo_sce, "pca.corrected")
  X_pca_batch <- X_pca
  
  for (b in names(cellids_sample)){
    batch_effect <- rnorm(ncol(X_pca), mean=0, sd = norm_sd)
    X_pca_batch[cellids_sample[[b]],] <- t(apply(X_pca_batch[cellids_sample[[b]],], 1, function(x) x + batch_effect))
  }
  
  reducedDim(embryo_sce, "pca_batch") <- X_pca_batch
  embryo_sce  
}

### METHODS ###

## Milo
# include function for running makenhoods
makeNhoods2 <- function(x, prop=0.1, k=21, d=30, refined=TRUE, reduced_dims="PCA", refinement_scheme = "reduced_dim") {
  if(is(x, "Milo")){
    message("Checking valid object")
    # check that a graph has been built
    if(!.valid_graph(miloR::graph(x))){
      stop("Not a valid Milo object - graph is missing. Please run buildGraph() first.")
    }
    graph <- miloR::graph(x)
    
    if(isTRUE(refined) & refinement_scheme == "reduced_dim"){
      X_reduced_dims  <- reducedDim(x, reduced_dims)
      if (d > ncol(X_reduced_dims)) {
        warning("Specified d is higher than the total number of dimensions in reducedDim(x, reduced_dims).
                        Falling back to using",ncol(X_reduced_dims),"dimensions\n")
        d <- ncol(X_reduced_dims)
      }
      X_reduced_dims  <- X_reduced_dims[,seq_len(d)]
      mat_cols <- ncol(x)
      match.ids <- all(rownames(X_reduced_dims) == colnames(x))
      if(!match.ids){
        stop("Rownames of reduced dimensions do not match cell IDs")
      }
    }
  } else if(is(x, "igraph")){
    
    if(isTRUE(refined) & refinement_scheme == "reduced_dim" & !is.matrix(reduced_dims)) {
      stop("No reduced dimensions matrix provided - required for refined sampling with refinement_scheme = reduced_dim.")
    }
    
    graph <- x
    
    if(isTRUE(refined) & refinement_scheme == "reduced_dim"){
      X_reduced_dims  <- reduced_dims
      mat_cols <- nrow(X_reduced_dims)
      if(is.null(rownames(X_reduced_dims))){
        stop("Reduced dim rownames are missing - required to assign cell IDs to neighbourhoods")
      }
    }
    
    if(isTRUE(refined) & refinement_scheme == "graph" & is.matrix(reduced_dims)){
      warning("Ignoring reduced dimensions matrix because refinement_scheme = graph was selected.")
    }
    
  } else{
    stop("Data format: ", class(x), " not recognised. Should be Milo or igraph.")
  }
  
  random_vertices <- .sample_vertices(graph, prop, return.vertices = TRUE)
  
  if (isFALSE(refined)) {
    sampled_vertices <- random_vertices
  } else if (isTRUE(refined)) {
    if(refinement_scheme == "reduced_dim"){
      sampled_vertices <- .refined_sampling(random_vertices, X_reduced_dims, k)
    } else if (refinement_scheme == "graph") {
      sampled_vertices <- .graph_refined_sampling(random_vertices, graph)
    } else {
      stop("When refined == TRUE, refinement_scheme must be one of \"reduced_dim\" or \"graph\".")
    }
  } else {
    stop("refined must be TRUE or FALSE")
  }
  
  sampled_vertices <- unique(sampled_vertices)
  
  if(is(x, "Milo")){
    nh_mat <- Matrix(data = 0, nrow=ncol(x), ncol=length(sampled_vertices), sparse = TRUE)
  } else if(is(x, "igraph")){
    nh_mat <- Matrix(data = 0, nrow=length(V(x)), ncol=length(sampled_vertices), sparse = TRUE)
  }
  # Is there an alternative to using a for loop to populate the sparseMatrix here?
  # if vertex names are set (as can happen with graphs from 3rd party tools), then set rownames of nh_mat
  v.class <- V(graph)$name
  
  if(is(x, "Milo")){
    rownames(nh_mat) <- colnames(x)
  } else if(is(x, "igraph")){
    if(is.null(v.class) & refinement_scheme == "reduced_dim"){
      rownames(nh_mat) <- rownames(X_reduced_dims)
    } else if(!is.null(v.class)){
      rownames(nh_mat) <- V(graph)$name
    }
  }
  
  for (X in seq_len(length(sampled_vertices))){
    nh_mat[unlist(neighborhood(graph, order = 1, nodes = sampled_vertices[X])), X] <- 1 #changed to include index cells
  }
  
  # need to add the index cells.
  colnames(nh_mat) <- as.character(sampled_vertices)
  if(is(x, "Milo")){
    nhoodIndex(x) <- as(sampled_vertices, "list")
    nhoods(x) <- nh_mat
    return(x)
  } else {
    return(nh_mat)
  }
}

.refined_sampling <- function(random_vertices, X_reduced_dims, k){
  message("Running refined sampling with reduced_dim")
  vertex.knn <-
    findKNN(
      X = X_reduced_dims,
      k = k,
      subset = as.vector(random_vertices),
      get.index = TRUE,
      get.distance = FALSE
    )
  
  nh_reduced_dims <- t(apply(vertex.knn$index, 1, function(x) colMedians(X_reduced_dims[x,])))
  
  # this function fails if rownames are not set
  if(is.null(rownames(X_reduced_dims))){
    warning("Rownames not set on reducedDims - setting to row indices")
    rownames(X_reduced_dims) <- as.character(seq_len(nrow(X_reduced_dims)))
  }
  
  colnames(nh_reduced_dims) <- colnames(X_reduced_dims)
  rownames(nh_reduced_dims) <- paste0('nh_', seq_len(nrow(nh_reduced_dims)))
  
  ## Search nearest cell to average profile
  # I have to do this trick because as far as I know there is no fast function to
  # search for NN between 2 distinct sets of points (here I'd like to search NNs of
  # nh_reduced_dims points among X_reduced_dims points). Suggestions are welcome
  all_reduced_dims <- rbind(nh_reduced_dims, X_reduced_dims)
  nn_mat <- findKNN(all_reduced_dims,
                    k = nrow(nh_reduced_dims) + 1,
                    subset = rownames(nh_reduced_dims))[["index"]]
  ## Look for first NN that is not another nhood
  nh_ixs <- seq_len(nrow(nh_reduced_dims))
  i = 1
  sampled_vertices <- rep(0, nrow(nn_mat))
  while (any(sampled_vertices <= max(nh_ixs))) {
    update_ix <- which(sampled_vertices <= max(nh_ixs))
    sampled_vertices[update_ix] <- nn_mat[update_ix, i]
    i <- i + 1
  }
  ## Reset indexes
  sampled_vertices <- sampled_vertices - max(nh_ixs)
  return(sampled_vertices)
}

.valid_graph <- function(x){
  # check for a valid graph
  if(isTRUE(is_igraph(x))){
    TRUE
  } else{
    FALSE
  }
}

.sample_vertices <- function(graph, prop, return.vertices=FALSE){
  # define a set of vertices and neihbourhood centers - extract the neihbourhoods of these cells
  random.vertices <- sample(V(graph), size=floor(prop*length(V(graph))))
  if(isTRUE(return.vertices)){
    return(random.vertices)
  } else{
    message("Finding neighbours of sampled vertices")
    vertex.list <- sapply(seq_len(length(random.vertices)), FUN=function(X) neighbors(graph, v=random.vertices[X]))
    return(list(random.vertices, vertex.list))
  }
}

.graph_refined_sampling <- function(random_vertices, graph){
  message("Running refined sampling with graph")
  random_vertices <- as.vector(random_vertices)
  X_graph <- set_vertex_attr(graph, "name", value = 1:length(V(graph)))
  refined_vertices <- lapply(seq_along(random_vertices), function(i){
    target_vertices <- unlist(neighborhood(X_graph, order = 1, nodes = random_vertices[i])) #get neighborhood of random vertex
    target_vertices <- target_vertices[-1] #remove first entry which is the random vertex itself
    rv_induced_subgraph <- induced_subgraph(graph = X_graph, vids = target_vertices)
    triangles <- count_triangles(rv_induced_subgraph)
    max_triangles <- max(triangles)
    max_triangles_indices <- which(triangles == max_triangles)
    #note - take first max_ego_index in the next line of code
    resulting_vertices <- V(rv_induced_subgraph)[max_triangles_indices]$name[1]
    return(resulting_vertices)
  }) %>% unlist() %>% as.integer()
  return(refined_vertices)
}

### as I cannot install the development branch, I will include code for the new version of testNhoods and spatialFDR here
testNhoods2 <- function(x, design, design.df,
                       fdr.weighting=c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"),
                       min.mean=0, model.contrasts=NULL, robust=TRUE, reduced.dim="PCA",
                       norm.method=c("TMM", "RLE", "logMS")){
  if(is(design, "formula")){
    model <- model.matrix(design, data=design.df)
    rownames(model) <- rownames(design.df)
  } else if(is(design, "matrix")){
    model <- design
    if(nrow(model) != nrow(design.df)){
      stop("Design matrix and model matrix are not the same dimensionality")
    }
    
    if(any(rownames(model) != rownames(design.df))){
      warning("Design matrix and model matrix dimnames are not the same")
      # check if rownames are a subset of the design.df
      check.names <- any(rownames(model) %in% rownames(design.df))
      if(isTRUE(check.names)){
        rownames(model) <- rownames(design.df)
      } else{
        stop("Design matrix and model matrix rownames are not a subset")
      }
    }
  }
  
  if(!is(x, "Milo")){
    stop("Unrecognised input type - must be of class Milo")
  } else if(.check_empty(x, "nhoodCounts")){
    stop("Neighbourhood counts missing - please run countCells first")
  }
  
  if(!any(norm.method %in% c("TMM", "logMS", "RLE"))){
    stop("Normalisation method ", norm.method, " not recognised. Must be either TMM, RLE or logMS")
  }
  
  if(!reduced.dim %in% reducedDimNames(x)){
    stop(reduced.dim, " is not found in reducedDimNames. Avaiable options are ", paste(reducedDimNames(x), collapse=","))
  }
  
  subset.counts <- FALSE
  if(ncol(nhoodCounts(x)) != nrow(model)){
    # need to allow for design.df with a subset of samples only
    if(all(rownames(model) %in% colnames(nhoodCounts(x)))){
      message("Design matrix is a strict subset of the nhood counts")
      subset.counts <- TRUE
    } else{
      stop("Design matrix (", nrow(model), ") and nhood counts (",
           ncol(nhoodCounts(x)), ") are not the same dimension")
    }
  }
  
  # assume nhoodCounts and model are in the same order
  # cast as DGEList doesn't accept sparse matrices
  # what is the cost of cast a matrix that is already dense vs. testing it's class
  if(min.mean > 0){
    if(isTRUE(subset.counts)){
      keep.nh <- rowMeans(nhoodCounts(x)[, rownames(model)]) >= min.mean
    } else{
      keep.nh <- rowMeans(nhoodCounts(x)) >= min.mean
    }
  } else{
    if(isTRUE(subset.counts)){
      keep.nh <- rep(TRUE, nrow(nhoodCounts(x)[, rownames(model)]))
    }else{
      keep.nh <- rep(TRUE, nrow(nhoodCounts(x)))
    }
  }
  
  if(isTRUE(subset.counts)){
    keep.samps <- intersect(rownames(model), colnames(nhoodCounts(x)[keep.nh, ]))
  } else{
    keep.samps <- colnames(nhoodCounts(x)[keep.nh, ])
  }
  
  if(any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) != rownames(model)) & !any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) %in% rownames(model))){
    stop("Sample names in design matrix and nhood counts are not matched.
             Set appropriate rownames in design matrix.")
  } else if(any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) != rownames(model)) & any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) %in% rownames(model))){
    warning("Sample names in design matrix and nhood counts are not matched. Reordering")
    model <- model[colnames(nhoodCounts(x)[keep.nh, keep.samps]), ]
  }
  
  if(length(norm.method) > 1){
    message("Using TMM normalisation")
    dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                   lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    dge <- calcNormFactors(dge, method="TMM")
  } else if(norm.method %in% c("TMM")){
    message("Using TMM normalisation")
    dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                   lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    dge <- calcNormFactors(dge, method="TMM")
  } else if(norm.method %in% c("RLE")){
    message("Using RLE normalisation")
    dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                   lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    dge <- calcNormFactors(dge, method="RLE")
  }else if(norm.method %in% c("logMS")){
    message("Using logMS normalisation")
    dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                   lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
  }
  
  dge <- estimateDisp(dge, model)
  fit <- glmQLFit(dge, model, robust=robust)
  if(!is.null(model.contrasts)){
    mod.constrast <- makeContrasts(contrasts=model.contrasts, levels=model)
    res <- as.data.frame(topTags(glmQLFTest(fit, contrast=mod.constrast),
                                 sort.by='none', n=Inf))
  } else{
    n.coef <- ncol(model)
    res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
  }
  
  res$Nhood <- as.numeric(rownames(res))
  message("Performing spatial FDR correction with", fdr.weighting[1], " weighting")
  mod.spatialfdr <- graphSpatialFDR(x.nhoods=nhoods(x),
                                    graph=miloR::graph(x),
                                    weighting=fdr.weighting,
                                    k=x@.k,
                                    pvalues=res[order(res$Nhood), ]$PValue,
                                    indices=nhoodIndex(x),
                                    distances=nhoodDistances(x),
                                    reduced.dimensions=reducedDim(x, reduced.dim))
  
  res$SpatialFDR[order(res$Nhood)] <- mod.spatialfdr
  res
}

graphSpatialFDR <- function(x.nhoods, graph, pvalues, k=NULL, weighting='k-distance',
                            reduced.dimensions=NULL, distances=NULL, indices=NULL){
  
  # Discarding NA pvalues.
  haspval <- !is.na(pvalues)
  
  if (!all(haspval)) {
    coords <- coords[haspval, , drop=FALSE]
    pvalues <- pvalues[haspval]
  }
  
  if(weighting[1] == "none"){
    return(rep(NA_real_, length(pvalues)))
  }
  # if the weighting vector length > 1 then just take the first element as this is the default
  # is this a bit hacky?
  if(length(weighting) > 1){
    weighting <- weighting[1]
  }
  
  if(weighting == "neighbour-distance"){
    if(!is.null(reduced.dimensions)){
      t.connect <- sapply(colnames(x.nhoods)[haspval],
                          FUN=function(PG) {
                            x.pcs <- reduced.dimensions[x.nhoods[, PG] > 0, ]
                            x.euclid <- as.matrix(dist(x.pcs))
                            x.distdens <- mean(x.euclid[lower.tri(x.euclid, diag=FALSE)])
                            return(x.distdens)})
    } else if(is.list(distances) & all(unlist(lapply(distances, class)) %in% c("matrix"))){
      t.connect <- unlist(lapply(distances, FUN=function(NHD) mean(rowMeans(NHD))))
    } else{
      stop("A matrix of reduced dimensions is required to calculate distances")
    }
  } else if(weighting == "max"){
    # do we have a distance matrix for the vertex cell to it's kth NN?
    if(!is.null(distances) & !is.null(indices)){
      # use distances first as they are already computed
      # compute the distance to the kth nearest neighbour
      # this is just the most distant neighbour
      if(class(distances) %in% c("matrix")){
        # find the distance to the kth nearest neighbour within the distance matrix
        t.connect <- unlist(lapply(indices, FUN=function(X) max(distances[X, ])))
      } else if(class(distances) %in% c("list")){
        t.connect <- unlist(lapply(indices, FUN=function(X) max(distances[[as.character(X)]])))
      } else{
        stop("Neighbourhood distances must be either a matrix or a list of matrices")
      }
    }
  }else if(weighting == "k-distance"){
    if(is.null(k)){
      stop("K must be non-null to use k-distance. Please provide a valid integer value")
    }
    
    # do we have a distance matrix for the vertex cell to it's kth NN?
    if(!is.null(distances) & !is.null(indices)){
      # use distances first as they are already computed
      # compute the distance to the kth nearest neighbour
      # this is just the most distant neighbour
      if(class(distances) %in% c("matrix")){
        # find the distance to the kth nearest neighbour within the distance matrix
        t.connect <- unlist(lapply(indices, FUN=function(X) distances[X, ][order(distances[X, ], decreasing=FALSE)[k]]))
      } else if(class(distances) %in% c("list")){
        # check if row names are set
        # distances is a list, so need to loop over slots to check for rownames
        null.names <- any(unlist(lapply(distances, FUN=function(RX) is.null(rownames(RX)))))
        if(isFALSE(null.names)){
          # nhood indices are _always_ numeric, distance slot names are strings - use the rownames of reducedDim to get the ID
          # this will need to be fixed properly across the whole code base
          t.dists <- lapply(indices, FUN=function(X) as.numeric((distances[[as.character(X)]])[rownames(reduced.dimensions)[X], ]))
          t.connect <- unlist(lapply(t.dists, FUN=function(Q) ((Q[Q>0])[order(Q[Q>0], decreasing=FALSE)])[k]))
        } else {
          # if row names are not set, extract numeric indices
          non.zero.nhoods <- which(x.nhoods != 0, arr.ind = TRUE)
          t.dists <- lapply(indices,
                            FUN=function(X) distances[[as.character(X)]][which(non.zero.nhoods[non.zero.nhoods[, 2] == which(indices == X),][, 1] == X),])
          t.connect <- unlist(lapply(t.dists, FUN=function(Q) (Q[Q>0])[order(Q[Q>0], decreasing=FALSE)[k]]))
        }
      } else{
        stop("Neighbourhood distances must be either a matrix or a list of matrices")
      }
    } else if(!is.null(reduced.dimensions) & !is.null(indices)){
      # find the kth NN and distance
      t.connect <- unlist(lapply(indices,
                                 FUN=function(X) max(findKNN(reduced.dimensions,
                                                             get.distance=TRUE,
                                                             subset=X, k=k)[["distance"]])))
    } else if(is.null(indices)){
      stop("No neighbourhood indices found - required to compute k-distance weighting")
    } else{
      stop("k-distance weighting requires either a distance matrix or reduced dimensions.")
    }
    
  } else if(weighting == "graph-overlap"){
    # no distance matrix is required here
    # compute overlap between neighborhoods
    if (!is.null(x.nhoods)) {
      intersect_mat <- crossprod(x.nhoods)
      diag(intersect_mat) <- 0
      t.connect <- unname(rowSums(intersect_mat))
    } else{
      stop("No neighborhoods found - please run makeNhoods first")
    }
    
  } else{
    stop("Weighting option not recognised - must be either k-distance, neighbour-distance, max or graph-overlap")
  }
  
  # use 1/connectivity as the weighting for the weighted BH adjustment from Cydar
  w <- 1/unlist(t.connect)
  w[is.infinite(w)] <- 0
  
  # Computing a density-weighted q-value.
  o <- order(pvalues)
  pvalues <- pvalues[o]
  w <- w[o]
  
  adjp <- numeric(length(o))
  adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
  adjp <- pmin(adjp, 1)
  
  if (!all(haspval)) {
    refp <- rep(NA_real_, length(haspval))
    refp[haspval] <- adjp
    adjp <- refp
  }
  return(adjp)
}

run_milo_red <- function(sce, condition_col, sample_col, reduced.dim="pca.corrected",
                         k=15, d=30, prop=0.1, returnMilo = TRUE,
                         batch_col=NULL, refined = TRUE, refinement_scheme = "reduced_dim"){
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Build graph neighbourhoods
  milo <- Milo(sce)
  milo <- buildGraph(milo, k=k, d=d, reduced.dim = reduced.dim)
  milo <- makeNhoods2(milo, prop = prop, k=k, d=d, reduced_dims = reduced.dim, 
                      refined = TRUE, refinement_scheme = "reduced_dim")
  ## Test DA
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
  milo <- calcNhoodDistance(milo, d=d, reduced.dim = reduced.dim)
  DA_results <- testNhoods2(milo, design = design, design.df = design_df, reduced.dim= reduced.dim)
  if (isTRUE(returnMilo)) {
    return(list(Milo=milo, DAres=DA_results))
  } else {
    DA_results 
  }
}

run_milo_tri <- function(sce, condition_col, sample_col, reduced.dim="pca.corrected",
                         k=15, d=30, prop=0.1, returnMilo = TRUE,
                         batch_col=NULL, refined = TRUE, refinement_scheme = "graph"){
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Build graph neighbourhoods
  milo <- Milo(sce)
  milo <- buildGraph(milo, k=k, d=d, reduced.dim = reduced.dim)
  milo <- makeNhoods2(milo, prop = prop, k=k, d=d, reduced_dims = reduced.dim, 
                      refined = TRUE, refinement_scheme = "graph", directed = FALSE, order = 2)
  ## Test DA
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
  #milo <- calcNhoodDistance(milo, d=d, reduced.dim = reduced.dim)
  DA_results <- testNhoods2(milo, design = design, fdr.weighting = "graph-overlap", design.df = design_df, reduced.dim= reduced.dim)
  if (isTRUE(returnMilo)) {
    return(list(Milo=milo, DAres=DA_results))
  } else {
    DA_results 
  }
}

milo2output <- function(milo, da_res, out_type="continuous", alpha=0.1){
  if (out_type=="continuous") { 
    da.cell.mat <- milo@nhoods %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.nhoods <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.nhoods.mat <- sapply(unique(da.nhoods), function(x) as.numeric(da.nhoods==x))
    da.cell.mat <- milo@nhoods %*% da.nhoods.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell
}

### RUN BENCHMARK ON SYNTHETIC LABELS ###

runDA <- function(sce, coldata, X_pca, 
                  method, 
                  condition_col='synth_labels', 
                  sample_col="synth_samples",
                  params = list(milo_red = list(k=50),
                                milo_tri = list(k=50),
                  ), d=30, out_type="label"
){
  ## Check that method name is in params
  if (!method %in% names(params)) {
    stop(paste("Specify parameters for method", method))
  }
  ## Check valid method
  if (!method %in% c("milo_red", "milo_tri")) {
    stop("Unrecognized method")
  }
  
  print("running DA")  
  ## Add reduced dim + coldata to sce
  colData(sce) <- DataFrame(coldata)
  reducedDim(sce, "pca_batch") <- as.matrix(X_pca)
  
  ## Run method
  ## Run milo
  if (method=="milo_red") {
    print("running milo_red")
    milo_res_red <- run_milo_red(sce, condition_col=condition_col, sample_col=sample_col,
                                 reduced.dim = "pca_batch", d=d, k=params$milo_red$k)
    out <- milo2output(milo_res_red$Milo, milo_res_red$DAres, out_type = out_type)
  } else if (method == "milo_tri"){
    print("running milo_tri")
    milo_res_tri <- run_milo_tri(sce, condition_col=condition_col, sample_col=sample_col,
                                 reduced.dim = "pca_batch", d=d, k=params$milo_tri$k)
    out <- milo2output(milo_res_tri$Milo, milo_res_tri$DAres, out_type = out_type)
  }
  
  ## Save results
  bm <- data.frame(bm_out=out)
  if (length(out) < ncol(sce)) {
    return(out)
  }
  bm$true_prob <- sce$Condition2_prob 
  bm$true <- sce$true_labels
  if (!is.null(sce$true_DA_clust)) {
    bm$true_clust <- sce$true_DA_clust
  }
  long_bm <- pivot_longer(bm, cols = bm_out, names_to='method', values_to="pred")
  long_bm[["method"]] <- method
  return(long_bm)
}

calculate_outcome <- function(long_bm){
  long_bm <- long_bm %>%
    mutate(outcome=case_when(true==pred & pred!="NotDA" ~ 'TP',
                             true!=pred & pred!="NotDA" ~ 'FP',
                             true!=pred & pred=="NotDA" ~ 'FN',
                             true==pred & pred=="NotDA"  ~ "TN"
    )) %>%
    group_by(method, outcome) %>%
    summarise(n=n()) %>%
    pivot_wider(id_cols=method, names_from=outcome, values_from=n, values_fill=0) 
  
  check_cols <- c("TP","FP","FN","TN") %in% colnames(long_bm)
  if (any(!check_cols)) {
    add_cols <- c("TP","FP","FN","TN")[!check_cols]
    for (col in add_cols) {
      long_bm[[col]] <- rep(0, nrow(long_bm))
    }
  }
  
  long_bm %>%
    mutate(TPR=TP/(TP+FN), FPR=FP/(FP+TN), TNR=TN/(TN+FP), FNR = FN/(FN+TP),
           FDR = FP/(TP+FP),
           Precision = TP/(TP+FP),
           Power = 1 - FNR,
           Accuracy = (TP + TN)/(TP + TN + FP + FN)
    )
}

## --- old functions --- ##

benchmark_da <- function(sce, condition_col='synth_labels', 
                         sample_col="synth_samples",
                         red_dim="pca.corrected",
                         params = list(milo = list(k=15),
                                       meld = list(k=15),
                                       daseq = list(k.vec=c(10,20,30, 40)),
                                       louvain = list(k=15),
                                       cydar = list(tol=5, downsample=3)
                         ),
                         d=30, out_type = "continuous"){
  ## Run milo
  milo_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                       reduced.dim = red_dim, d=d, k=params$milo$k)
  milo_out <- milo2output(milo_res$Milo, milo_res$DAres, out_type = out_type)
  ## Run milo controlling for batch
  milo_batch_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                             reduced.dim = red_dim, d=d, k=params$milo$k, batch_col = "synth_batches")
  milo_batch_out <- milo2output(milo_batch_res$Milo, milo_batch_res$DAres, out_type = out_type)
  
  ## Run DAseq
  daseq_res <- run_daseq(sce, k.vec=params$daseq$k.vec, condition_col, 
                         reduced.dim = red_dim, d=d)
  daseq_out <- daseq2output(sce, daseq_res, out_type = out_type)
  ## Run MELD
  meld_res <- run_meld_reticulate(sce, condition_col=condition_col, sample_col=sample_col,
                                  reduced.dim = red_dim, d=d, k=params$meld$k)
  meld_out <- meld2output(meld_res, out_type = out_type)
  ## Run louvain
  louvain_res <- run_louvain(sce, condition_col=condition_col, sample_col=sample_col,
                             reduced.dim = red_dim, d=d, k=params$louvain$k)
  louvain_out <- louvain2output(louvain_res, out_type = out_type)
  ## Collect results + true labels
  bm <- data.frame(milo=milo_out, milo_batch=milo_batch_out,
                   daseq=daseq_out, meld=meld_out, louvain=louvain_out)
  bm$true_prob <- sce$Condition2_prob 
  bm$true <- sce$true_labels
  if (!is.null(sce$true_DA_clust)) {
    bm$true_clust <- sce$true_DA_clust
  }
  long_bm <- pivot_longer(bm, cols = c(milo, milo_batch, daseq, meld, louvain), names_to='method', values_to="pred") 
  return(long_bm)
}


# Given a data_embedding, sample a simplex to weight each dimension to get a
# new PDF for each condition
# a: scaling coefficient in logit (lower a --> less extreme probabilities) 
.make_pdf <- function(data_embedding, a=0.2){
  # Create an array of values that sums to 1
  n_components = ncol(data_embedding)
  data_simplex = sort(runif(n = n_components-1))
  data_simplex = c(0, data_simplex, 1)
  data_simplex = diff(data_simplex)
  data_simplex = sample(data_simplex)
  # Weight each embedding component by the simplex weights
  sort_axis = rowSums(data_embedding * data_simplex)
  
  # Pass the weighted components through a logit
  pdf = 1/(1+exp(- a * sort_axis))
  if (sample(c(TRUE, FALSE), 1)){ pdf = 1 - pdf }
  return(pdf)  
}

## Smooth probabilities over KNN graph (to avoid having clusters with opposite sign DA)
.knn_smoothing <- function(graph, cond_probability, k=15, redDim='pca.corrected', d=10){
  # X_red_dim = reducedDim(sce, redDim)[,1:d]
  # graph = buildKNNGraph(t(X_red_dim), k = k)  
  adj = get.adjacency(graph)
  smooth_cond_probability <- (adj %*% cond_probability)/rowSums(adj)
  smooth_cond_probability
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels <- function(sce, # SingleCellExperiment obj
                                 knn_graph, # for knn smoothing of probability values
                                 redDim='pca.corrected', # embedding to use to simulate differential abundance
                                 n_conditions=2, # number of conditions to simulate
                                 n_components=10, # number of components of embedding to use
                                 n_replicates=3, # number of replicates per condition
                                 n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                 seed=42){
  data_embedding = reducedDim(sce, redDim)[,1:n_components]
  set.seed(seed)
  # embedding data must be mean-centered
  data_embedding = t(scale(t(data_embedding), scale=FALSE))
  
  # Randomly flip sign of each embedding dimension
  data_embedding = apply(data_embedding, 2, function(x)  x*sample(c(-1, 1), 1) )
  
  conditions = paste0("Condition", 1:n_conditions)
  cond_probability = sapply(1:(length(conditions)-1), function(x) .make_pdf(data_embedding))
  
  # KNN Smoothing to avoid regions of the graph with opposite labels
  cond_probability = .knn_smoothing(knn_graph, cond_probability, redDim=redDim, d=n_components)
  
  # Normalize to sum to 1 for each cell
  # cond_probability <- t(apply(cond_probability, 1, function(x) x/sum(abs(x))))
  cond_probability = cbind(cond_probability, 1 - rowSums(cond_probability))
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  synth_samples <- paste0(synth_labels, "_", replicates)
  names(batches) <- sort(unique(synth_samples))
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}