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

million_cells <- readRDS("/nfs/research/marioni/alsu/hubmap_metaRef/data/lung/sce_mapped.Rds")
print("million cells imported")

million_cells <- million_cells[,apply(reducedDim(million_cells, "pca.corrected"), 1, function(x) !all(is.na(x)))]
#million_cells <- runUMAP(million_cells, dimred = "pca.corrected", name = 'umap')

# paste makenhood devel function 
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

# insert devel testnhoods
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

## These are utility functions not meant to be exposed to the user

.check_empty <- function(x, attribute){
    # check if a Milo object slot is empty or not
    x.slot <- slot(x, attribute)

    if(is.list(x.slot) & names(slot(x, "graph")) == "graph"){
        return(length(x.slot[[1]]) > 0)
    } else if(is.list(x.slot) & is.null(names(x.slot))){
        return(length(x.slot))
    } else if(any(class(x.slot) %in% c("dgCMatrix", "dsCMatrix", "ddiMatrix", "matrix"))){
        return(sum(rowSums(x.slot)) == 0)
    }
}

.check_binary <- function(x){
    # check if a matrix is binary or not
    sum.zeros <- sum(x == 0)
    sum.ones <- sum(x == 1)
    n.comps <- nrow(x) * ncol(x)

    return(sum(c(sum.zeros, sum.ones)) == n.comps)
}

.neighborsToKNNGraph <- function(nn, directed=FALSE) {
    start <- as.vector(row(nn))
    end <- as.vector(nn)
    interleaved <- as.vector(rbind(start, end))

    if (directed) {
        g <- make_graph(interleaved, directed=TRUE)

    } else {
        g <- make_graph(interleaved, directed=FALSE)
        g <- simplify(g, edge.attr.comb = "first")
    }
    g
}

.set_reduced_dims <- function(x, value, slot.x=NULL, rdim=NULL){
    x <- updateObject(x)
    content <- slot(x, slot.x)

    if(slot.x == "nhoodReducedDim"){

        if(!is.null(rdim)){
            content[[rdim]] <- value
            x@nhoodReducedDim <- content
        } else{
            stop("No reduced dimensionality slot provided")
        }
    }else{
        stop("replacement method not implemented for ", slot)
    }

    x
}

#### nhood adjacency matrix function

.build_nhood_adjacency <- function(nhoods, overlap=1){
    nh_intersect_mat <- Matrix::crossprod(nhoods)
    nh_intersect_mat[nh_intersect_mat < overlap] <- 0

    rownames(nh_intersect_mat) <- colnames(nhoods)
    colnames(nh_intersect_mat) <- colnames(nhoods)
    return(nh_intersect_mat)
}

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
                                million_cells3 <- makeNhoods2(million_cells2, prop = 0.1, k = k, d=30, refined = TRUE, reduced_dims = "pca.corrected", refinement_scheme = "graph");
                                million_cells4 <- countCells(million_cells3, meta.data = as.data.frame(colData(million_cells3)), samples="sample2");
                                da_results <- testNhoods2(million_cells4, reduced.dim= "pca.corrected", design = ~ health_status, design.df = design, fdr.weighting = "graph-overlap")})
    
    k_val[j] <- k
    graph[j] <- a[3]
    j <- j + 1
}

time.df <- cbind(k_val, graph)
time.df <- data.frame(time.df)
write_csv(time.df, "/nfs/research/marioni/akluzer/speed_benchmarks/time_millioncells.csv")
