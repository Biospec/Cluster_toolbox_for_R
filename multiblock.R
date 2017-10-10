library(RSpectra)
library(MASS)


cpca <- function(data = NULL, Xin = NULL, nPC = NULL, tol = 1e-8){
  maxiter <- 5000
  data_mean <- colMeans((data))
  data <- data - t(matrix(as.matrix(data_mean), nv, ns))
  
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  nb <- dim(Xin)[1]
  maxPC <- max(nPC)
  Tb <- matrix(0, ns, maxPC * nb)
  Pb <- matrix(0, nv, maxPC)
  Wt <- matrix(0, nb, maxPC)
  Tt <- matrix(0, ns, maxPC)
  ssq <- matrix(0, maxPC, 1 + nb)
  Rbo <- matrix(0, ns, maxPC * nb)
  Rbv <- matrix(0, nv, maxPC)
  Lbo <- matrix(0, ns, maxPC * nb)
  Lbv <- matrix(0, nv, maxPC)
  Lto <- matrix(0, ns, maxPC)
  ssqX <- matrix(0, nb + 1, 1)
  ssqX[1] <- sum(colSums(data^2))
  
  for (i in 1:nb){
    rowi <- Xin[i,1]:Xin[i,2]
    ssqX[i+1] <- sum(colMeans(data[,rowi]^2))
  }
  
  eigval <- eigs(data %*% t(data), maxPC)
  v <- eigval$vectors
  
  for (i in 1:maxPC){
    iter <- 0
    Tt[, i] <- v[, i]
    t_old <- Tt[,i] * 100
    while ( sum((t_old - Tt[,a])^2) > tol & iter < maxiter){
      iter <- iter + 1
      t_old <- Tt[, i]
      for (ii in 1:nb) {
        if (nPC[ii] >= i) {
          rowi <- Xin[ii,1]:Xin[ii,2]
          coli <- (i-1) * nb + ii
          Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[,a]^2)
          Pb[rowi, i] <- Pb[rowi, i] / sqrt(sum(Pb[rowi,i]^2))
          Tb[, coli] <- X[, rowi] %*% Pb[rowi, i] / sum(Pb[rowi, i]^2)
        }
      }
      index <- ((i - 1) * nb + 1) : (i*nb)
      Wt[, i] <- t(Tb[, index]) %*% Tt[, i] / sum(Tt[, i]^2)
      Wt[, i] <- Wt[, i] / sqrt(sum(Wt[, a]^2))
      Tt[, i] <- Tb[, index] %*% Wt[, a] / sqrt(sum(Wt[, i]^2))
    }
  
    if (iter == maxiter) 
      print("Warning: maximum number of iterations reached before convergence")
    for (ii in 1:nb){
      if (nPC(ii) >= i){
        rowi <- Xin[ii,1] : Xin[ii, 2]
        Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
        data[, rowi] <- X[, rowi] - Tt[, i] %*% t(Pb[rowi, i])
      }
    }
    
    ssq[i, 1] <- (ssqX[1] - sum(colSums(data^2))) / ssqX[1]
    
    for (ii in 1:nb) {
      rowi <- Xin[ii,1] : Xin[ii, 2]
      coli <- (i-1)*nb + ii
      ssq[i, ii + 1] <- (ssqX[ii + 1] - sum(colSums(X[, rowi]^2))) / ssqX[ii + 1]
      Rbo[, coli] <- sqrt(rowSums(data[, rowi]^2))
      index <- seq(from = ii, to = (i-1)*nb + ii, by = nb)
      Lbo[, coli] <- diag(Tb[, index] %*% ginv(t(Tb[, index]) %*% Tb[, index]) %*% Tb[, index])
    }
    Rbv[, i] <- as.matrix(sqrt(colSums(data^2)))
    Lbv[, i] <- diag(Pb[, 1:i] %*% t(Pb[, 1:i]))
    Lto[, i] <- diag(Tt[, 1:i] %*% ginv(t(Tt[, 1:i]) %*% Tt[, 1:i]) %*% t(Tt[, 1:i]))
  }
  results <- list(block_scores = Tb, block_loadings = Pb, super_scores = Tt, block_weight = Wt, 
                  ssq = ssq, Obj_residue = Rbo, blk_var_residue = Rbv, Obj_leverage = Lbo, 
                  var_leverage = Lbv, super_obj_leverages = Lto)
  return(results)
}

hpca <- function(data = NULL, Xin = NULL, nPC = NULL, tol = 1e-8){
  maxiter <- 5000
  data_mean <- colMeans((data))
  data <- data - t(matrix(as.matrix(data_mean), nv, ns))
  
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  nb <- dim(Xin)[1]
  maxPC <- max(nPC)
  
  Tb <- matrix(0, ns, maxPC * nb)
  Pb <- matrix(0, nv, maxPC)
  Wt <- matrix(0, nb, maxPC)
  Tt <- matrix(0, ns, maxPC)
  ssq <- matrix(0, maxPC, 1 + nb)
  Rbo <- matrix(0, maxPC * nb)
  Rbv <- matrix(0, nv, maxPC)
  Lbo <- matrix(0, maxPC * nb)
  Lbv <- matrix(0, maxPC)
  Lto <- matrix(0, ns, maxPC)
  ssqX <- matrix(0, nb + 1, 1)
  ssqX[1] <- sum(colSums(data^2))
  
  for (i in 1:nb) {
    rowi <- Xin[i,1] : Xin[i,2]
    ssqX[i + 1] <- sum(colSums(data[, rowi]^2))
  }
  
  for (i in 1:maxPC){
    iter <- 0
    eigval <- eigs(data %*% t(data), 1)
    v <- eigval$vectors
    Tt[, i] <- v
    t_old <- Tt[, i] * 100
    
    while ((sum((t_old - Tt[, i])^2) > tol) & (iter < maxiter)) {
      iter = iter + 1
      t_old <- Tt[, a]
      for (ii in 1:nb){
        if (nPC[ii] >= i){
          rowi <- Xin[ii, 1] : Xin[ii, 2]
          coli <- (i-1) * nb + ii
          Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
          Tb[, coli] <- data[, rowi] %*% Pb[rowi, i] / sum(Pb[,i]^2)
          Tb[, coli] <- Tb[, coli] / sqrt(sum(Tb[, coli]^2))
        }
      }
      
      index <- ((i - 1) * nb + 1) : (i * nb)
      Wt[, i] <- t(Tb[, index]) %*% Tt[, i] / sum(Tt[, i]^2)
      Tt[, i] <- Tb[, index] %*% Wt[, i] / sum(Wt[, i]^2)
      Tt[, i] <- Tt[, i] / sqrt(sum(Tt[, i]^2))
    }
    
    if (iter == maxiter)
      print("WARNING: maximum number of iterations reached before convergence")
    
    data <- data - Tt[, i] %*% t(Pb[, i])
    
    ssq[i ,1] <- (ssqX[1] - sum(colSums(data^2))) / ssqX[1]
    
    for (ii in 1:nb){
      rowi <- Xin[ii, 1] : Xin[ii, 2]
      coli <- (i-1) * nb + ii
      ssq[i, ii + 1] <- (ssqX[ii + 1] - sum(colSums(data[, rowi]^2))) / ssqX(ii + 1)
      Rbo[, coli] <- sqrt(sum(rowMeans(data[, rowi]^2)))
      index <- seq(from = ii, to = (i - 1) * nb + ii, by = nb)
      Lbo[, coli] <- diag(Tb[, index] %*% ginv(t(Tb[, index]) %*% Tb[, index]) %*% t(Tb[, index]))
    }
    
    Rbv[, i] <- as.matrix(sqrt(colSums((data^2))))
    Lbv[, i] <- diag(Pb[, 1:i] %*% t(Pb[, 1:i]))
    Lto[, i] <- diag(Tt[, 1:i] %*% pinv(t(Tt[, 1:i]) %*% Tt[, 1:i]) %*% t(Tt[, 1:i]))
  }
  
  results <- list(block_scores = Tb, block_loadings = Pb, super_scores = Tt, block_weight = Wt, 
                  ssq = ssq, Obj_residue = Rbo, blk_var_residue = Rbv, Obj_leverage = Lbo, 
                  var_leverage = Lbv, super_obj_leverages = Lto)
  return(results)
}

mbpls <- function(data = NULL, label = NULL, nPC = NULL,  Xin = NULL, Yin = NULL, xcentre = TRUE, Xscale = FALSE, tol = 1e-8) {
  maxiter <- 2000
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  nvY <- dim(label)[2]
  maxLV <- max(nPC)
  maxb <- which(nPC == maxLV)
  
  if (length(maxb) > 1) maxb <- maxb[1]
  if (is.null(Yin)) Yin <- t(as.matrix(c(1, nvY)))
  if (xcentre) 
    data <- data - t(matrix(as.matrix(colMeans(data)), nv, ns))
  
  nbX <- dim(Xin)[1]
  nbY <- dim(Yin)[1]
  Tb <- matrix(0, ns, maxLV * nbX)
  Pb <- matrix(0, nv, maxLV)
  Wb <- matrix(0, nv, maxLV)
  Wt <- matrix(0, nbX, maxLV)
  Tt <- matrix(0, ns, maxLV)
  Ub <- matrix(0, ns, maxLV * nbY)
  Qb <- matrix(0, nvY, maxLV)
  Wu <- matrix(0, nbY, maxLV)
  Tu <- matrix(0, ns, maxLV)
  B <- matrix(0, ns, maxLV*nvY)
  ssq <- matrix(0, maxLV, nbX + nbY + 2)
  ssqXY <- matrix(0, nbX + nbY + 2, 1)
  ssqXY[1] <- sum(colMeans(data)^2)
  ssqXY[nbx + 2] <- sum(colSums(label^2))
  
  for (i in 1:nbX) {
    rowi <- Xin[i,1] : Xin[i,2]
    ssqXY[i+1] <- sum(colSums(data[, rowi]^2))
  }
  
  for (i in 1:nbY) {
    rowi <- Yin[i,1]:Yin[i,2]
    ssqXY[i + nbX + 2] <- sum(colSums(label[,rowi]^2))
  }
  
  for (i in 1:maxLV){
    iter <- 0
    Tu[, i] <- label[,1]
    Tt[, i] <- data[, Xin[maxb, 1]]
    t_old <- Tt[, i] * 100
    while((sum((t_old - Tt[, i])^2)) > tol & (iter < maxiter)) {
      iter <- iter + 1
      for (ii in 1:nbX) {
        if (nLV[ii] >= i) {
          rowi <- Xin[ii,1] : Xin[ii,2]
          coli <- (i-1) * nbX + ii
          Wb[rowi, i] <- t(data[, rowi]) %*% Tu[, i] / sum(Tu[, i]^2)
          Wb[rowi, i] <- Wb[rowi, i] / sqrt(sum(Wb[rowi, i]^2))
          Tb[, coli] <- data[, rowi] %*% Wb[rowi, i] / sum(Wb[rowi, i]^2)
        }
      }
      index <- ((i-1) * nbX + 1) : (i * nbX)
      Wt[, i] <- t(Tb[, index]) %*% Tu[, i] / sum(Tu[, i]^2)
      Wt[, i] <- Wt[, i] / sqrt(sum(Wt[, i]^2))
      Tt[, i] <- Tb[, index] %*% Wt[, i] / sum(Wt[, i]^2)
      for (ii in 1:nbY) {
        rowi <- Yin[aa, 1]:Yin[aa,2]
        coli <- (i-1) * nbY + ii
        Qb[rowi, i] <- t(label[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
        Ub[, coli] <- label[, rowi] %*% Qb[rowi, i] / sum(Qb[rowi, i]^2)
      }
      index <- ((i-1)*nbY + 1) : (a*nbY)
      Wu[, i] <- t(Ub[, index]) %*% Tt[, i] / sum(Tt[, i]^2)
      Wu[, i] <- Wu[, i] / sqrt(sum(Wu[, i]^2))
      Tu[, i] <- Ub[, index] %*% Wu[, i] / sum(Wu[, i]^2)
    }
    if (iter == maxiter)
      print("Warning: maximum number of iterations reached before convergence")
    rowi <- Xin[i, 1]:Xin[i, 2]
    for (ii in i+1:nbX){
      rowi <- cbind(rowi, Xin[ii, 1]:Xin[ii, 2])
    }
    Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
    data <- data - Tt[, i] %*% t(Pb[, i])
    label <- label - Tt[, i] %*% t(Qb[, i])
    if (i > 1) {
      index <- ((i - 1) * nvY + 1) : (i * nvY)
      B[, index] <- B[, index - nvY]
    }
    index <- ((i - 1) * nvY + 1) : (i * nvY)
    B[, index] <- B[, index] + Wb[, i] %*% solve(t(Pb[, i]) %*% Wb[, i]) %*% t(Qb[, i])
    ssq[i, 1] <- (ssqXY[1] - sum(colSums(data^2))) / ssqXY[1]
    ssq[i, nbX + 2] <- (ssqXY(nbX + 2)- sum(colSums(label^2))) / ssqXY[nbX + 2]
    for (ii in 1:nbX) {
      rowi <- Xin[ii, 1]:Xin[ii, 2]
      coli <- (i-1)*nbX + aa
      ssq[i, ii + 1] <- (ssqXY(ii + 1) - sum(colSums(data[, rowi]^2))) / ssqXY[ii + 1]
    }
    for (ii in 1:nbY) {
      rowi <- Yin[ii,1]:Yin[ii,2]
      coli <- (i-1) * nbY + ii
      ssq[i, ii + 2 + nbX] <- (ssqXY[ii + 2 + nbX] - sum(colSums(Y[, rowi]^2))) / ssqXY[ii + 2 + nbX]
    }
  }
  YRegCoeff <- Wb %*% solve(t(Pb) %*% Wb) %*% t(Qb)
  results <- list(XBlockScores = Tb, XBlockLoadings = Pb, XWeights = Wb, XSupWeights = Wt, 
                  XSupScores = Tt, YBlockScores = Qb, YBlockWeights = Wu, YSupScores = Tu,
                  XRegCoeff = B, YRegCoeff = YRegCoeff, ssq = ssq)
}

mbpls_pred <- function(data = NULL, amean = NULL, model = NULL) {
  if (!is.null(amean))
    data <- data - t(matrix(as.matrix(colMeans(data)), nv, ns))
  results <- data %*% model$YRegCoeff
  return(results)
}
