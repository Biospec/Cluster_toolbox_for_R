library(pracma)

pcdfa <- function(x, label = NULL, no_pc = 10, maxfac = NULL) {
    if (is.null(label))
        label <- row.names(x)
    unique_label <- unique(label)
    no_class <- length(unique_label)
    if (no_class < 2) stop("Number of classes must > 1")
    if (maxfac > no_class - 1 || is.null(maxfac)) {
        print("No. of DF is missing or too many DFs to be extracted, set it to group - 1")
        maxfac <- no_class - 1
    }

    xmean <- colMeans(x)
    pca_results <- prcomp(x)
    x <- pca_results$x[,1:no_pc]
    pc_loadings <- pca_results$rotation[,1:no_pc]

    mx <- colMeans(x)

    K <- x[label == unique_label[1],]
    ns <- dim(K)[1]
    nv <- dim(K)[2]
    if (ns > 1) zz <- colMeans(K) else zz <- K
    A <- K - t(matrix(mx, nv, ns))
    C <- K - t(matrix(zz, nv, ns))
    Tt <- t(A) %*% A
    W <- t(C) %*% C
    for (i in 2:no_class) {
        K = x[label == unique_label[i],]
        ns <- dim(K)[1]
        nv <- dim(K)[2]
        if (ns > 1) zz <- colMeans(K) else zz <- K
        A <- K - t(matrix(mx, nv, ns))
        C <- K - t(matrix(zz, nv, ns))
        Tt <- Tt + t(A) %*% A
        W <- W + t(C) %*% C
    }
    B <- Tt - W
    invW <- solve(W)
    P <- invW %*% B
    P_eigen <- eigen(P)
    P_eigen$values <- Re(P_eigen$values[1:maxfac])
    P_eigen$vectors <- Re(P_eigen$vectors[, 1:maxfac])
    V <- P_eigen$vectors
    U <- x %*% V
    loadings <- pc_loadings %*% V
    results <- list(scores = U, loadings = loadings, eigenvalues = P_eigen$values, xmean = xmean)
}

pcdfa_pred <- function(x, model = NULL) {
    x <- t(apply(x, 1, "-", model$xmean))
    prediction <- x %*% model$loadings
    return(prediction)
}


pls <- function(X = NULL, Y = NULL, lv = NULL, Xmc = TRUE, Ymc = TRUE, Xscale = FALSE, Yscale = FALSE) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    ns <- dim(X)[1]
    nv <- dim(X)[2]
    nvY <- dim(Y)[2]
    
    if (is.null(lv)) lv <- min(c(ns, nv))
    
    if (Xmc) {
      Xmean <- colMeans(X)
      X <- X - t(matrix(Xmean, nv, ns))
    } else
      Xmean <- NULL
    
    if (Ymc) {
      if (nvY > 1){
        Ymean <- colMeans(Y)
        Y <- Y- t(matrix(Ymean, nvY, ns))
      } else {
        Ymean <- mean(Y)
        Y <- Y - Ymean
      }
    } else
      Ymean <- NULL
    
    if (is.function(Xscale)) {
      Xscale <- apply(X,2, scale[1])
      X <- X / t(matrix(Xscale, nv, ns))
    }else if (Xscale) {
      Xscale <- apply(X, 2, sd)
      X <- X / t(matrix(Xscale, nv, ns))
    } else
      Xscale <- NULL
    
    if (is.function(Yscale)) {
      Yscale <- apply(Y, 2, scale[2])
    }else if (Yscale) {
      Yscale <- apply(Y, 2, sd)
    } else
      Yscale <- NULL
    
    b <- matrix(0, lv, 1)
    Tt <- matrix(0, ns, lv)
    P <- matrix(0, nv, lv)
    W <- matrix(0, nv, lv)
    Q <- matrix(0, nvY, lv)
    for (i in 1:lv) {
        XY <- t(X) %*% Y
        usv <- svd(XY)
        q <- usv$v[, 1]
        w <- usv$u[, 1]
        ts <- X %*% w
        p <- t(ts) %*% X / as.double((t(ts) %*% ts))
        p <- t(p)
        u <- Y %*% q
        b[i] <- 1 / (t(ts) %*% ts) * t(ts) %*% u
        X <- X - ts %*% t(p)
        Y <- Y - b[i] * ts %*% t(q)
        Tt[, i] <- ts
        W[, i] <- w
        Q[, i] <- q
        P[, i] <- p
    }
    if (lv == 1)
        R <- Q * as.double(b)
    else    
        R <- Q %*% diag(as.vector(b))
    B <- W %*% solve(t(P) %*% W) %*% t(R)
    results <- list(Xscores = Tt, Xloadings = P, Weights = W, Yloadings = Q, b = b, B = B, 
                    Xmean = Xmean, Xscale = Xscale, Ymean = Ymean, Yscale = Yscale)
    return(results)
}

plspred <- function(model = NULL, X = NULL) {
    Xmean <- model$Xmean
    Xscale <- model$Xscale
    Ymean <- model$Ymean
    Yscale <- model$Yscale
    
    ns <- dim(X)[1]
    nv <- dim(X)[2]
    nvY <- dim(model$Yloadings)[1]
    B <- model$B
    
    if (!is.null(Xmean)) X = X - t(matrix(Xmean, nv, ns))
    if (!is.null(Xscale)) X = X / t(matrix(Xscale, nv, ns))
    pred <- X %*% B

    if (!is.null(Ymean)) pred <- pred + t(matrix(Ymean, nvY, ns))
    if (!is.null(Yscale)) pred <- pred * t(matrix(Yscale, nvY, ns))
    results <- pred
    return(results)
}

plsda <- function(X = NULL, Y = NULL, lv = NULL, ...) {
    Y <- as.matrix(Y)
    ns <- dim(Y)[1]
    nvY <- dim(Y)[2]
    if (nvY == 1) {
        unique_Y <- unique(Y)
        no_class <- length(unique_Y)
        Y_new <- matrix(0, ns, no_class)
        for (i in 1:no_class) {
            Y_new[Y == unique_Y[i], i] <- 1
        }
    } else {
        Y_new <- Y
    }
    results <- pls(X, Y_new, lv, ...)
    return(results)
}

plsdapred <- function(model = NULL, X = NULL) {
    raw_pred <- plspred(model, X)
    ns <- dim(raw_pred)[1]
    nYv <- dim(raw_pred)[2]
    pred <- matrix(0, ns, 1)    
    for (i in 1:ns) {
        pred[i] <- which(raw_pred[i,] == max(raw_pred[i,]))
    }
    results <- pred
    return(results)
}





