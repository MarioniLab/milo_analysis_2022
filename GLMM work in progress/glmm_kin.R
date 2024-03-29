#' Perform differential abundance testing using a NB-generalised linear mixed model
#' 
#' This function will perform DA testing on all nhoods using a negative binomial generalised linear mixed model
#'
#' @param X A matrix containing the fixed effects of the model.
#' @param Z A matrix containing the random effects of the model.
#' @param y A matrix containing the observed phenotype over each neighborhood. 
#' @param REML A logical value denoting whether REML (Restricted Maximum Likelihood) should be run. Default is TRUE.
#' @param random.levels A list describing the random effects of the model, and for each, the different unique levels.
#' @param glmm.control A list containing parameter values specifying the theta tolerance of the model and the maximum number of iterations to be run.
#' @param dispersion A scalar value for the dispersion of the negative binomial.
#' 
#' @details 
#' This function runs a negative binomial generalised linear mixed effects model. If mixed effects are detected in testNhoods, 
#' this function is run to solve the model.
#'
#' @importMethodsFrom Matrix %*%
#' @importFrom stats runif
#' @importFrom Matrix Matrix solve
#' @export
runGLMMkin <- function(X, Z, y, random.levels=NULL, REML=TRUE,
                       glmm.control=list(theta.tol=1e-6, max.iter=50),
                       dispersion = NULL, Kin=NULL){
    
    print("running Alice's Kin in R")
    
    # model components
    # X - fixed effects model matrix
    # Z - random effects model matrix
    # y - observed phenotype
    theta.conv <- glmm.control[["theta.tol"]] # convergence for the parameters
    max.hit <- glmm.control[["max.iter"]]
    r <- dispersion
    
    # OLS for the betas is usually a good starting point for NR                    
    curr_beta <- solve((t(X) %*% X)) %*% t(X) %*% log(y + 1)  
    rownames(curr_beta) <- colnames(X)
    
    # create full Z with expanded random effect levels
    full.Z <- initializeFullZ.R(Z)
    
    # random value initiation from runif
    curr_u <- matrix(runif(colnames(full.Z)), ncol=1)
    rownames(curr_u) <- colnames(full.Z)
    
    # compute sample variances of the us
    curr_sigma <- Matrix(runif(ncol(Z), 0, 1), ncol=1, sparse = TRUE)
    #rownames(curr_sigma) <- colnames(Z)
    rownames(curr_sigma) <- "Genetic" #####here
    
    random.levels <- c(list("Genetic"=seq_len(nrow(curr_u)))) #####here
    
    #compute variance-covariance matrix G
    curr_G <- initialiseGR(cluster_levels=random.levels, sigmas=curr_sigma, Kin=Kin) #####here
    
    #try adding Z* = ZL #####here
    #L <- t(chol(curr_G))
    #full.Z <- full.Z %*% L 
    
    #create a single variable for the thetas
    curr_theta <- do.call(rbind, list(curr_beta, curr_u))
    
    #compute mu.vec using inverse link function                   
    mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
    
    theta_diff <- rep(Inf, nrow(curr_theta))
    sigma_diff <- Inf
    
    #compute variance-covariance matrix G
    G_inv <- solve(curr_G)
    
    conv.list <- list()
    iters <- 1
    meet.conditions <- !((all(theta_diff < theta.conv)) & (sigma_diff < theta.conv) | iters >= max.hit)
    
    while(meet.conditions){
        
        print(curr_sigma)
        print(curr_beta)
        #print(curr_G[1:5, 1:5])
        
        #compute all matrices - information about them found within their respective functions
        D <- computeD(mu=mu.vec)
        D_inv <- solve(D)
        y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u, y=y)
        V <- computeV(mu=mu.vec, r=r)
        W <- computeW(D_inv=D_inv, V=V)
        W_inv <- solve(W)
        V_star <- computeV_star(full.Z=full.Z, curr_G=curr_G, W=W)
        V_star_inv <- solve(V_star)
        #V_partial <- computeV_partial(full.Z=full.Z, random.levels=random.levels, curr_sigma=curr_sigma)
        V_partial <- computeV_partial_kin(Kin=Kin, full.Z=full.Z, random.levels=random.levels) #####here
        #V_partial <- computeV_partial_stable(full.Z=full.Z, random.levels=random.levels, curr_sigma=curr_sigma)
        
        #---- First estimate variance components with Newton Raphson procedure ---#
        if (isFALSE(REML)) {
            score_sigma <- sigmaScore(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, random.levels=random.levels)
            information_sigma <- sigmaInformation(V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels)
        } else if (isTRUE(REML)) {
            P <- computeP_REML(V_star_inv=V_star_inv, X=X)
            PV <- computePV(V_partial=V_partial, P=P)
            loglihood <- computeloglihood(n=nrow(X), V_star=V_star, y_star=y_star, X=X, curr_beta=curr_beta, V_star_inv=V_star_inv)
            score_sigma <- sigmaScoreREML(PV=PV, V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, P=P, random.levels=random.levels)
            information_sigma <- sigmaInformationREML(PV=PV, V_star_inv=V_star_inv, V_partial=V_partial, P=P, random.levels=random.levels)
        }
        sigma_update <- FisherScore(score_vec=score_sigma, hess_mat=information_sigma, theta_hat=curr_sigma, random.levels=random.levels)
        sigma_diff <- abs(sigma_update - curr_sigma)
        
        # update sigma, G, and G_inv
        curr_sigma <- sigma_update
        #rownames(curr_sigma) <- colnames(Z)
        rownames(curr_sigma) <- "Genetic" #####here
        curr_G <- initialiseGR(cluster_levels=random.levels, sigmas=curr_sigma, Kin=Kin) #####here
        G_inv <- solve(curr_G)
        
        #---- Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u---####
        theta_update <- solve_equations(X=X, W_inv=W_inv, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star) 
        #theta_update <- solve_equations_kin3(X=X, curr_sigma=curr_sigma, W_inv=W_inv, Kin=Kin, curr_G=curr_G, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star)
        
        if (isTRUE(theta_update)) {
            stop("Hessian is computationally singular - cannot solve GLMM")
        }
        
        theta_diff <- abs(theta_update - curr_theta)
        
        # update B, u and mu_vec to determine new values of score and hessian matrices
        curr_theta <- theta_update
        rownames(curr_theta) <- c(colnames(X), colnames(full.Z))
        curr_beta <- curr_theta[colnames(X), , drop=FALSE]
        curr_u <- curr_theta[colnames(full.Z), , drop=FALSE]
        mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
        
        if (any(is.infinite(mu.vec))) {
            stop("Estimates increasing to infinity - cannot solve GLMM.")
        }
        
        conv.list[[iters]] <- list("FE"=as.vector(curr_beta),
                                   "RE"=as.vector(curr_u),
                                   "Sigma"=as.vector(curr_sigma),
                                   "Theta.Converged"=theta_diff < theta.conv,
                                   "Sigma.Converged"=sigma_diff < theta.conv,
                                   "theta_diff"=theta_diff,
                                   "sigma_diff"=sigma_diff,
                                   "Iters"=iters,
                                   "loglihood"=loglihood)
        
        iters <- iters + 1
        meet.conditions <- !((all(theta_diff < theta.conv)) & (all((sigma_diff) < theta.conv))| iters >= max.hit) 
    }
    
    SE <- calculateSE(X=X, full.Z=full.Z, W_inv=W_inv, G_inv=G_inv)
    Zscore <- calculateZscore(curr_beta=curr_beta, SE=SE)
    df <- Satterthwaite_df.R(X=X, PV=PV, SE=SE, REML=REML, W_inv=W_inv, full.Z=full.Z, curr_sigma=curr_sigma, curr_beta=curr_beta, random.levels=random.levels, V_partial=V_partial, V_star_inv=V_star_inv, G_inv=G_inv)
    Pvalue <- computePvalue(Zscore=Zscore, df=df)
    
    converged <- ((all(theta_diff < theta.conv)) & (all(abs(sigma_diff) < theta.conv)))
    final.list <- list("FE"=as.vector(curr_beta),
                       "RE"=as.vector(curr_u),
                       "Sigma"=as.vector(curr_sigma),
                       "Theta.Converged"=theta_diff < theta.conv,
                       "Sigma.Converged"=sigma_diff < theta.conv,
                       "converged"=converged,
                       "Iters"=iters,
                       "Dispersion"=r, 
                       "SE"=SE,
                       "Zscore"=Zscore,
                       "df" = df,
                       "pvalue"=Pvalue,
                       "ystar" = y_star,
                       "conv.list"=conv.list)
    
    return(final.list)
}

computeloglihood <- function(n, V_star, y_star, X, curr_beta, V_star_inv){
    loglihood <- -(n/2)*log(2*pi)-0.5*det(V_star)-0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% (y_star - X %*% curr_beta)
}

#' @importMethodsFrom Matrix %*%
#' @export
computeW <- function(D_inv=D_inv, V=V){
    W = D_inv %*% V %*% D_inv
    return(W)
}

#' @importFrom Matrix Matrix
#' @export
computeV <- function(mu, r){
    # compute diagonal matrix of variances
    v.vec <- ((mu**2/r)) + mu
    V <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(V) <- v.vec
    return(V)
}

#' @importFrom Matrix Matrix
#' @export
computeD <- function(mu=mu.vec){
    # D is diag(mu_i)
    D <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(D) <- mu
    return(D)
}

#' @importMethodsFrom Matrix %*%
#' @export
computeV_star <- function(full.Z=full.Z, curr_G=curr_G, W=W){
    V_star = full.Z %*% curr_G %*% t(full.Z) + W
    return(V_star)
}

#' @importMethodsFrom Matrix %*%
#' @export
computey_star <- function(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u, y=y){
    y_star <- ((X %*% curr_beta) + (full.Z %*% curr_u)) + D_inv %*% (y - exp((X %*% curr_beta) + (full.Z %*% curr_u)))
    return(y_star)
}

#' @importMethodsFrom Matrix %*%
#' @export
computeV_partial <- function(full.Z=full.Z, random.levels=random.levels, curr_sigma=curr_sigma){
    V_partial_vec <- list()
    indx <- c(0, cumsum(lengths(random.levels)))
    for (i in 1:length(random.levels)) {
        Z.temp <- full.Z[ , (indx[i]+1):(indx[i+1])]
        V_partial_vec[[i]] <- Z.temp %*% t(Z.temp)
    }
    return(V_partial_vec)
}

computeV_partial_kin <- function(Kin=Kin, full.Z=full.Z, random.levels=random.levels){
    V_partial_vec <- list()
    indx <- c(0, cumsum(lengths(random.levels)))
    for (i in 1:length(random.levels)) {
        Z.temp <- full.Z[ , (indx[i]+1):(indx[i+1])]
        V_partial_vec[[i]] <- Z.temp %*% Kin %*% t(Z.temp)
    }
    return(V_partial_vec)
}

computeV_partial_stable <- function(full.Z=full.Z, curr_sigma=curr_sigma, random.levels=random.levels){
    V_partial_vec <- list()
    indx <- c(0, cumsum(lengths(random.levels)))
    for (i in 1:length(random.levels)) {
        Z.temp <- full.Z[ , (indx[i]+1):(indx[i+1])]
        V_partial_vec[[i]] <- 0.5 * (1/sqrt(curr_sigma[1])) * Z.temp
    }
    return(V_partial_vec)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
sigmaScore <- function(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, random.levels=random.levels){
    score_vec <- Matrix(0, nrow = length(V_partial), ncol = 1, sparse = TRUE)
    for (i in 1:length(V_partial)){
        score_vec[i,1] <- -0.5*matrix.trace(V_star_inv %*% V_partial[[i]]) + 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% V_partial[[i]] %*% V_star_inv %*% (y_star - X %*% curr_beta)
    }
    return(score_vec)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
sigmaInformation <- function(V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels) {
    info_vec <- Matrix(sapply(1:length(random.levels), function(i){
        0.5*matrix.trace(V_star_inv %*% V_partial[[i]] %*% V_star_inv %*% V_partial[[i]])}), ncol =1, sparse = TRUE)
    return(info_vec)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
computePV <- function(V_partial=V_partial, P=P){
    PV <- list()
    for (i in 1:length(V_partial)) {
        PV[[i]] <- P %*% V_partial[[i]]
    }
    return(PV)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
sigmaScoreREML <- function(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, P=P, random.levels=random.levels, PV=PV){
    score_vec <- Matrix(0, nrow = length(V_partial), ncol = 1, sparse = TRUE)
    for (i in 1:length(V_partial)) {
        score_vec[i,1] <- -0.5*matrix.trace(PV[[i]]) + 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% V_partial[[i]] %*% V_star_inv %*% (y_star - X %*% curr_beta)
    }
    return(score_vec)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
sigmaInformationREML <- function(V_star_inv=V_star_inv, V_partial=V_partial, P=P, random.levels=random.levels, PV=PV) {
    info_vec <- Matrix(sapply(1:length(random.levels), function(i){
        0.5*matrix.trace(PV[[i]] %*% PV[[i]])}), ncol=1, sparse = TRUE)
    return(info_vec)
}

#' @importMethodsFrom Matrix %*%
#' @export
computeP_REML <- function(V_star_inv=V_star_inv, X=X) {
    P <- V_star_inv - V_star_inv %*% X %*% solve(t(X) %*% V_star_inv %*% X) %*% t(X) %*% V_star_inv
    return(P)
}

#' @export
FisherScore <- function(score_vec, hess_mat, theta_hat, random.levels, lambda=1e-5, det.tol=1e-10, cond.tol=1e-15){
    # sequentially update the parameter using the Newton-Raphson algorithm
    # theta ~= theta_hat + hess^-1 * score
    # this needs to be in a direction of descent towards a minimum
    theta_new <- theta_hat + 1/(hess_mat) * score_vec
    return(theta_new)
}

#' @importFrom Matrix sparseMatrix diag
#' @export
initialiseGR <- function(cluster_levels, sigmas, Kin=NULL){
    # construct the correct size of G given the random effects and variance components
    # names of cluster_levels and columns of Z must match
    # the independent sigmas go on the diagonal and the off-diagonal are the crossed/interactions
    # sigmas must be named
    sum.levels <- sum(unlist(lapply(cluster_levels, length)))
    G <- sparseMatrix(i=sum.levels, j=sum.levels, repr="C", x=0L)
    dimnames(G) <- list(unlist(cluster_levels), unlist(cluster_levels))
    i <- j <- 1
    m <- ncol(G)
    
    for(x in seq_len(nrow(sigmas))){
        x.q <- length(cluster_levels[[rownames(sigmas)[x]]])
        if(is.null(Kin)){
            diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] # is this sufficient to transform the sigma to the model scale?
        } else{
            if(rownames(sigmas[x, , drop=FALSE]) %in% c("Genetic")){
                print("using Kin")
                G[c((m-x.q+1):m), c((m-x.q+1):m)] <- Kin * sigmas[x, ]
                #diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] * Kin
            }else{
                diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] # is this sufficient to transform the sigma to the model scale?
            }
        }
        
        i <- j <- i+x.q
    }
    return(as.matrix(G))
}

#' @importFrom Matrix Matrix
#' @export
initializeFullZ.R <- function(Z) {
    full.Z <- matrix(,nrow=nrow(Z), ncol = 0)
    for (i in 1:ncol(Z)) {
        temp.Z <- Matrix(table(seq_along(1:nrow(Z)), Z[,i]), sparse = TRUE)
        full.Z <- cbind(full.Z, temp.Z)
    }
    return(full.Z)
}

#' @importFrom Matrix Matrix
#' @export
solve_equations <- function(X=X, W_inv=W_inv, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star){
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv
    
    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(full.Z) %*% W_inv %*% y_star))
    theta_update <- tryCatch({solve(LHS) %*% RHS}
                             , error = function(e){ 
                                 exit <- TRUE
                             }
    )
    return(theta_update)
}

solve_equations_kin <- function(X=X, W_inv=W_inv, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star){
    message("using PCG")
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv
    
    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(full.Z) %*% W_inv %*% y_star))
    
    theta_update <- tryCatch({solve(LHS) %*% RHS}
                             , error = function(e){
                                 exit <- TRUE
                             })
    
    #theta_update <- cgsolve(A = as.matrix(LHS), b = as.matrix(RHS))
    #theta_update <- pcgsolve(A = as.matrix(LHS), b = as.matrix(RHS), preconditioner = "Jacobi")
    
    
    return(theta_update)
}

solve_equations_kin2 <- function(X=X, W_inv=W_inv, curr_sigma=curr_sigma, full.Z=full.Z, curr_G=curr_G, G_inv=G_inv, Kin=Kin, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star){
    message("using kin2 cholesky")
    
    L <- t(chol(Kin))
    S <- full.Z %*% L
    
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% S
    LowerLeft <- t(S) %*% W_inv %*% X
    LowerRight <- t(S) %*% W_inv %*% S + (diag(1, nrow = nrow(Kin), ncol = ncol(Kin)) %*% solve(diag(curr_sigma[1,], ncol(Kin), nrow(Kin))))
    
    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(S) %*% W_inv %*% y_star))
    
    # theta_update <- tryCatch({solve(LHS) %*% RHS}
    #                           , error = function(e){ 
    #                               exit <- TRUE
    #                           })
    
    #theta_update <- cgsolve(A = as.matrix(LHS), b = as.matrix(RHS))
    theta_update <- pcgsolve(A = as.matrix(LHS), b = as.matrix(RHS), preconditioner = "Jacobi")
    theta_update[3:nrow(theta_update)] <- L %*% theta_update[3:nrow(theta_update)]
    
    return(theta_update)
}

solve_equations_kin3 <- function(X=X, W_inv=W_inv, curr_sigma=curr_sigma, full.Z=full.Z, curr_G=curr_G, G_inv=G_inv, Kin=Kin, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star){
    message("using kin3 eigenvalue decomposition")
    
    L <- t(chol(Kin))
    S <- full.Z %*% L
    D_inv <- solve(diag(eigen(Kin)$values))
    
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% S
    LowerLeft <- t(S) %*% W_inv %*% X
    LowerRight <- t(S) %*% W_inv %*% S + (D_inv %*% solve(diag(curr_sigma[1,], ncol(Kin), nrow(Kin))))
    
    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(S) %*% W_inv %*% y_star))
    
    # theta_update <- tryCatch({solve(LHS) %*% RHS}
    #                          , error = function(e){ 
    #                              exit <- TRUE
    #                          })
    
    #theta_update <- cgsolve(A = as.matrix(LHS), b = as.matrix(RHS))
    theta_update <- pcgsolve(A = as.matrix(LHS), b = as.matrix(RHS), preconditioner = "Jacobi")
    theta_update[3:nrow(theta_update)] <- L %*% theta_update[3:nrow(theta_update)]
    
    return(theta_update)
}

#' @export
matrix.trace <- function(x){
    # check is square matrix first
    x.dims <- dim(x)
    if(x.dims[1] != x.dims[2]){
        stop("matrix is not square")
    } else{
        return(sum(diag(x)))
    }
}

conj_grad <- function(A, b, x){
    r <- b - A %*% x
    p <- r
    rsold <- t(r) %*% r
    
    for (i in 1:length(b)){
        Ap <- A %*% p 
        alpha = rsold / (t(p) %*% Ap)
        x <- x + alpha[1,1] * p 
        r <- r - alpha[1,1] * Ap
        rsnew <- t(r) %*% r
        if (sqrt(rsnew) < 1e-10) {
            return(x)
        }
        p <- r + (rsnew / rsold)[1,1] * p 
        rsold <- rsnew
    }
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix solve diag
#' @export
calculateSE <- function(X=X, full.Z=full.Z, W_inv=W_inv, G_inv=G_inv) {
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv
    vcov <- solve(UpperLeft - UpperRight %*% solve(LowerRight) %*% LowerLeft)
    se <- sqrt(diag(vcov))
    return(se)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
calculateZscore <- function(curr_beta=curr_beta, SE=SE) {
    Zscore <- as.matrix(curr_beta)/SE
    return(Zscore)
}

#' @importFrom stats pt
#' @export
computePvalue <- function(Zscore=Zscore, df=df) {
    pval <- 2*pt(abs(Zscore), df, lower.tail=FALSE)
    return(pval)
}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix solve diag
#' @importFrom numDeriv jacobian
#' @export
Satterthwaite_df.R <- function(X=X, PV=PV, SE=SE, REML=REML, W_inv=W_inv, full.Z=full.Z, curr_sigma=curr_sigma, curr_beta=curr_beta, random.levels=random.levels, V_partial=V_partial, V_star_inv=V_star_inv, G_inv=G_inv) {
    
    ###---- first calculate g = derivative of C with respect to sigma ----
    function_jac <- function(x, X.fun=as.matrix(X), W_inv.fun=as.matrix(W_inv), full.Z.fun=as.matrix(full.Z)) {
        UpperLeft <- t(X.fun) %*% W_inv.fun %*% X.fun
        UpperRight <- t(X.fun) %*% W_inv.fun %*% full.Z.fun
        LowerLeft <- t(full.Z.fun) %*% W_inv.fun %*% X.fun
        LowerRight <- t(full.Z.fun) %*% W_inv.fun %*% full.Z.fun
        n <- length(random.levels)
        diag(LowerRight) <- diag(LowerRight) + rep(1/x, times=lengths(random.levels)) #when extending to random slopes, this needs to be changed to a matrix and added to LowerRight directly
        C <- solve(UpperLeft - UpperRight %*% solve(LowerRight) %*% LowerLeft)
    }
    
    jac <- numDeriv::jacobian(func=function_jac, x=as.vector(curr_sigma))
    jac_list <- lapply(1:ncol(jac), function(i)
        array(jac[, i], dim=rep(length(curr_beta), 2))) #when extending to random slopes, this would have to be reformatted into list, where each element belongs to one random effect
    
    #next, calculate V_a, the asymptotic covariance matrix of the estimated covariance parameters
    #given by formula below
    P <- computeP_REML(V_star_inv=V_star_inv, X=X)
    V_a <- matrix(NA, nrow=length(random.levels), ncol=length(random.levels))
    for (i in 1:length(random.levels)) {
        for (j in 1:length(random.levels)) {
            V_a[i,j] <- 2*(1/(matrix.trace(PV[[i]] %*% PV[[j]])))
        }
    }
    
    df <- rep(NA, length(curr_beta))
    for (i in 1:length(curr_beta)) {
        jac_var_beta <- unlist(lapply(lapply(jac_list, diag), `[[`, i))
        denom <- t(jac_var_beta) %*% (V_a) %*% jac_var_beta #g' Va g
        df[i] <- 2*((SE[i]^2)^2)/denom
    }
    return(as.matrix(df))
}