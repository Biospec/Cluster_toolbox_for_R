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

plsboots <- function(data = NULL, label = NULL, rep_idx = NULL, iter = 1000, lv = NULL, type = "r", perm = FALSE) {
    if (is.null(rep_idx)) rep_idx <- 1:dim(data)[1]
    switch(type,
        r = plsr_boots(data, label, rep_idx, iter, lv, perm),
        d = plsda_boots(data, label, rep_idx, iter, lv), perms)
}

plsr_boots <- function(data, label, rep_idx, iter, lv, perm) {
    if (is.null(lv)) {
        lv <- min(dim(data))
    }
    data <- as.matrix(data)
    label <- as.matrix(label)
    rep_idx <- as.matrix(rep_idx)
    nvY <- dim(label)[2]

    opt_rmsecv <- matrix(0, nvY, iter)
    opt_Q2 <- matrix(0, nvY, iter)

    rmsep <- matrix(0, iter, nvY)
    Q2p <- matrix(0, iter, nvY)
    

    predictions <- vector("list", iter)
    if (perm) {
        rmsep_perm <- matrix(0, iter, nvY)
        Q2p_perm <- matrix(0, iter, nvY)    
        predictions_perm <- vector("list", iter)
    }
    for (i in 1:iter) {
        print(i)
        boot_idx <- boots(rep_idx)
        data_trn <- data[boot_idx$trn_idx,]
        label_trn <- as.matrix(label[boot_idx$trn_idx,])
        rep_trn <- as.matrix(rep_idx[boot_idx$trn_idx,])
        data_tst <- data[boot_idx$tst_idx,]
        ns <- dim(data_tst)[1]
        label_tst <- as.matrix(label[boot_idx$tst_idx,])
        Q2_cv <- matrix(0, nvY, lv)
        rmsecv <- matrix(0, nvY, lv)
        for (ii in 1:lv) {
            cv_results <- plsr_crossval(data = data_trn, label = label_trn, rep_idx = rep_trn, lv = ii)
            Q2_cv[,ii] <- cv_results$Q2
            rmsecv[,ii] <- cv_results$rmsecv
        }
        if (nvY == 1) {
            opt_LV <- which(Q2_cv == max(Q2_cv))
            if (length(opt_LV) > 1) opt_LV <- opt_LV[1]
            opt_rmsecv[i] <- rmsecv[opt_LV]
            opt_Q2[i] <- Q2_cv[opt_LV]
        } else {
            meanQ2 <- colMeans(Q2_cv)
            opt_LV <- which(meanQ2 == max(meanQ2))
            if (length(opt_LV) > 1) opt_LV <- opt_LV[1]
            opt_rmsecv[, i] <- rmsecv[, opt_LV]
            opt_Q2[,i] <- Q2_cv[, opt_LV]
        }
                
        amean <- t(as.matrix(colMeans(data_trn)))
        cmean <- t(as.matrix(colMeans(label_trn)))
        data_trn <- t(apply(data_trn, 1, "-", amean))
        if (nvY == 1)
            label_trn <- as.matrix(as.vector(label_trn) - as.double(cmean))
        else
            label_trn <- t(apply(label_trn, 1, "-", cmean))
        plsmodel <- pls(data_trn, label_trn, lv = opt_LV)
        predC <- plspred(data_tst, model = plsmodel, aMean = amean, cMean = cmean)
        
        if (nvY == 1) {
            rmsep[i,] <- sqrt(sum((label_tst - predC$pred) ^ 2) / ns)
            Q2p[i,] <- 1 - sum((label_tst - predC$pred) ^ 2) / sum((label_tst - mean(label_tst)) ^ 2)
            predictions[[i]] <- predC$pred
        }
        else {
            for (ii in 1:nvY) {
                rmsep[i, ii] <- sqrt(sum((label_tst[, ii] - predC$pred[, ii]) ^ 2) / ns)
                Q2p[i,ii] <- 1 - sum((label_tst[,ii] - predC$pred[,ii])^2) / sum((label_tst[,ii] - mean(label_tst[,ii]))^2)
            }
            predictions[[i]] <- predC$pred
        }
        if (perm) {
            perm_idx <- reps_perm(rep_trn)
            data_perm <- data_trn[perm_idx,]
            plsmodel_perm <- pls(data_perm, label_trn, lv = opt_LV)
            predC_perm <- plspred(data_tst, model = plsmodel_perm, aMean = amean, cMean = cmean)
            if (nvY == 1) {
                rmsep_perm[i,] <- sqrt(sum((label_tst - predC_perm$pred) ^ 2) / ns)
                Q2p_perm[i,] <- 1 - sum((label_tst - predC_perm$pred) ^ 2) / sum((label_tst - mean(label_tst)) ^ 2)
                predictions_perm[[i]] <- predC_perm$pred
            }
            else {
                for (ii in 1:nvY) {
                    rmsep_perm[i, ii] <- sqrt(sum((label_tst[, ii] - predC_perm$pred[, ii]) ^ 2) / ns)
                    Q2p_perm[i, ii] <- 1 - sum((label_tst[, ii] - predC_perm$pred[, ii]) ^ 2) / sum((label_tst[, ii] - mean(label_tst[, ii])) ^ 2)
                }
                predictions_perm[[i]] <- predC_perm$pred
            }
            
        }
    }
    if (perm) {
        if (dim(Q2p)[2] == 1)
            pval <- length(which(Q2p <= Q2p_perm)) / length(Q2p)
        else
            pval <- length(which(rowMeans(Q2p) <= rowMeans(Q2p_perm))) / dim(Q2p)[1]
        results <- list(predictions = predictions, predictions_perm = predictions_perm, rmsep = rmsep,
            rmsep_perm = rmsep_perm, Q2p = Q2p, Q2p_perm = Q2p_perm, Q2cv = opt_Q2, rmsecv = opt_rmsecv,
            pval = pval)
    }
    else
        results <- list(predictions = predictions, rmsep = rmsep, rmsecv = opt_rmsecv, Q2p = Q2p, Q2cv = opt_Q2)
    return(results)
}

plsda_boots <- function(data, label, rep_idx, iter=1000, lv, perm=FALSE) {
  if (is.null(lv)) {
    lv <- min(dim(data))
  }
  data <- as.matrix(data)
  label <- as.matrix(label)
  rep_idx <- as.matrix(rep_idx)
  ns <- dim(label)[1]
  nvY <- dim(label)[2]
  if (length(nvY) == 1) {
    unique_class = unique(label)
    nc <- length(unique_class)
    label_da <- matrix(0, ns, nc)
    for (i in 1:length(unique_class))
      label_da[label == unique_class[i], i] = 1
  } 
  else {
    label_da <- label
    label <- matrix(0, ns, 1)
    nc <- dim(label_da)[2]
    for (i in 1:ns)
      label[i] <- which(label_da[i,] == max(label_da[i,]))
  }
  
  ccr_p <- matrix(0, iter, 1)
  optccr_cv <- matrix(0, iter, 1)
  conf_mat <- replicate(iter, matrix(0,nc,nc))
  
  
  predictions <- vector("list", iter)
  
  if (perm) {
    ccr_perm <- matrix(0, iter, 1)
    conf_mat_perm <- replicate(iter, matrix(0, nc ,nc))
  }
  for (i in 1:iter) {
      print(i)        
      while (TRUE) {
            boot_idx <- boots(rep_idx)
            data_trn <- data[boot_idx$trn_idx,]
            label_trn <- label[boot_idx$trn_idx,]
            if (length(unique(label_trn)) == length(unique(label))) break
      }     
      rep_trn <- as.matrix(rep_idx[boot_idx$trn_idx,])
      data_tst <- data[boot_idx$tst_idx,]
      ns <- dim(data_tst)[1]
      label_tst <- label[boot_idx$tst_idx,]
      ccr_cv <- matrix(0, lv, 1)
    
    for (ii in 1:lv) {
      cv_results <- plsda_crossval(data = data_trn, label = label_trn, rep_idx = rep_trn, lv = ii)
      ccr_cv[ii] <- cv_results$ccr
    }
    
    opt_LV <- which(ccr_cv == max(ccr_cv))
    if (length(opt_LV) > 1) opt_LV <- opt_LV[1]

    amean <- t(as.matrix(colMeans(data_trn)))
    
    data_trn <- t(apply(data_trn, 1, "-", amean))   
    
    plsmodel <- plsda(data_trn, label_trn, lv = opt_LV)
    pred <- plsdapred(data_tst, model = plsmodel, label_known = label_tst, amean = amean)
    ccr_p[i] <- pred$ccr
    conf_mat[, , i] <- pred$conf_mat
    predictions[[i]] <- pred$predictions
    
    if (perm) {
      perm_idx <- reps_perm(rep_trn)
      data_perm <- data_trn[perm_idx,]
      plsmodel_perm <- plsda(data_perm, label_trn, lv = opt_LV)
      pred_perm <- plsdapred(data_tst, model = plsmodel_perm, label_known = label_tst, amean = amean)
      ccr_perm[i] <- pred_perm$ccr
      conf_mat_perm[, , i] <- pred_perm$conf_mat
    }
  }
  if (perm) {
    pval <- length(which(ccr_p <= ccr_perm))
    results <- list(predictions = predictions, ccr_p = ccr_p, ccr_perm = ccr_perm, 
                    ccr_cv = ccr_cv, pval = pval, conf_mat = conf_mat, 
                    conf_mat_perm = conf_mat_perm)
  }
  else
    results <- list(predictions = predictions, ccr_p = ccr_p, ccr_cv = ccr_cv, conf_mat = conf_mat)
  return(results)
}


reps_perm <- function(reps = NULL) {
    no_samples <- length(reps)
    unique_reps <- unique(reps)
    no_reps <- length(unique_reps)
    reps_after_perm <- unique_reps[randperm(no_reps)]
    perm_idx <- matrix(0, no_samples, 1)
    idx <- 1
    for (i in 1:no_reps) {
        reps_found <- which(reps == reps_after_perm[i])
        no_reps_found <- length(reps_found)
        perm_idx[idx:(idx + no_reps_found - 1)] <- reps_found
        idx <- idx + no_reps_found
    }
    return(perm_idx)
}

plsr_crossval <- function(data, label, rep_idx = NULL, lv = NULL, k = NULL) {
    label <- as.matrix(label)
    nvY <- dim(label)[2]
    rep_idx <- as.matrix(rep_idx)
    unique_rep <- unique(rep_idx)
    ns <- length(unique_rep)
    if (nvY > 1){
      rmsecv <- matrix(0, nvY, 1)
      Q2 <- matrix(0, nvY, 1)
    } 
    
    if (is.null(lv)) lv <- min(dim(data))
    
    if (is.null(k)) {
        if (length(unique_rep) < 14) k <- ns else k <- 7
        }    
    step <- round(ns / k)
    predCval <- matrix(0, dim(label)[1], dim(label)[2])
    for (i in seq(1, ns, by = step)) {
        rep_val <- unique_rep[i:min(c(i + step - 1, ns))]
        val_idx <- which(!is.na(match(rep_idx, rep_val)))
        data_val <- data[val_idx,]
        label_val <- as.matrix(label[val_idx,])
        data_trn <- data[-val_idx,]
        label_trn <- as.matrix(label[-val_idx,])
        amean <- t(as.matrix(colMeans(data_trn)))
        cmean <- t(as.matrix(colMeans(label_trn)))
        data_trn <- t(apply(data_trn, 1, "-", amean))
        if (nvY > 1)
            label_trn <- t(apply(label_trn, 1, "-", cmean))
        else
            label_trn <- as.matrix(apply(label_trn, 1, "-", cmean))
        plsmodel <- pls(data_trn, label_trn, lv = lv)
        predC <- plspred(data_val, plsmodel, aMean = amean, cMean = cmean)
        predCval[val_idx,] <- predC$pred
    }
    if (nvY == 1) {
        rmsecv <- sqrt(sum((label - predCval) ^ 2) / ns)
        Q2 <- 1 - sum((label - predCval) ^ 2) / sum((label - mean(label)) ^ 2)
    } else {
        for (ii in 1:nvY) {
            rmsecv[ii] <- sqrt(sum((label[, ii] - predCval[, ii]) ^ 2) / ns)
            Q2[ii] <- 1 - sum((label[, ii] - predCval[,ii]) ^ 2) / sum((label[, ii] - mean(label[, ii])) ^ 2)
        }
    }
      
    results <- list(predictions = predCval, rmsecv = rmsecv, Q2 = Q2)
    return(results)
}

plsda_crossval <- function(data, label, rep_idx = NULL, lv = NULL, k = NULL) {
    if (is.null(lv)) lv <- min(dim(data))
    label <- as.matrix(label)
    unique_rep <- unique(rep_idx)
    ns <- length(unique_rep)

    if (is.null(k)) {
        if (length(unique_rep) < 14) k <- ns else k <- 7
    }    
    step <- round(ns / k)    
    nsY <- dim(label)[1]
    nvY <- dim(label)[2]
    predLabel <- matrix(0, nsY, 1)

    if (nvY > 1) {
        label_da <- label
        label <- matrix(0, nsY, 1)
        nc <- dim(label_da)[2]
        label <- as.matrix(apply(matrix(as.logical(label_da), nsY, nvY), 1, "which"))
    } else {
        unique_class <- unique(label)
        nc <- length(unique_class)
        label_da <- matrix(0, nsY, nc)
        label_new <- matrix(0, nsY,1)
        for (i in 1:nc) {
            label_da[label == unique_class[i], i] <- 1
            label_new[label == unique_class[i],] <- i
        }
        label <- label_new
        rm(label_new)
    }
    for (i in seq(1, ns, by = step)) {
        rep_val <- unique_rep[i:min(c(i + step - 1, ns))]
        val_idx <- which(!is.na(match(rep_idx, rep_val)))
        data_val <- data[val_idx,]
        label_val <- as.matrix(label_da[val_idx,])
        data_trn <- data[-val_idx,]
        label_trn <- as.matrix(label_da[-val_idx,])
        amean <- t(as.matrix(colMeans(data_trn)))        
        data_trn <- t(apply(data_trn, 1, "-", amean))        
        plsmodel <- plsda(data_trn, label_trn, lv = lv)
        predC <- plsdapred(data_val, model = plsmodel, amean = amean)
        predLabel[val_idx,] <- predC$label_predict
    }

    ccr <- length(which(predLabel == label)) / length(label)
    conf_mat <- matrix(0, nc, nc)
    for (i in 1:nc) {
        for (ii in 1:nc) {
            conf_mat[i, ii] <- length(which(label == i & predLabel == ii)) / length(which(label == i))
        }
    }
    results <- list(ccr = ccr, conf_mat = conf_mat, predictions = predLabel)
}



pls <- function(X = NULL, Y = NULL, lv = NULL) {
    ns <- dim(X)[1]
    nv <- dim(X)[2]
    nvY <- dim(Y)[2]
    if (is.null(lv)) lv <- min(c(ns, nv))
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
    results <- list(Xscores = Tt, Xloadings = P, Weights = W, Yloadings = Q, b = b, B = B)
    return(results)
}

plspred <- function(data = NULL, model = NULL, aMean = NULL, cMean = NULL, ascale = NULL, cscale = NULL) {
    ns <- dim(data)[1]
    nv <- dim(data)[2]
    nvY <- dim(model$Yloadings)[1]
    B <- model$B
    
    if (!is.null(aMean)) data = data - t(matrix(aMean, nv, ns))
    if (!is.null(ascale)) data = data / t(matrix(ascale, nv, ns))
    pred <- data %*% B

    if (!is.null(cMean)) pred <- pred + t(matrix(cMean, nvY, ns))
    if (!is.null(cscale)) pred <- pred * matrix(cscale, nvY, ns)
    results <- list(pred = pred, B = B)
    return(results)
}

plsda <- function(data = NULL, label = NULL, lv = NULL) {
    label <- as.matrix(label)
    ns <- dim(label)[1]
    nvY <- dim(label)[2]
    if (nvY == 1) {
        unique_label <- unique(label)
        no_class <- length(unique_label)
        label_new <- matrix(0, ns, no_class)
        for (i in 1:no_class) {
            label_new[label == unique_label[i], i] <- 1
        }
    } else {
        label_new <- label
    }
    results <- pls(data, label_new, lv)
    return(results)
}

plsdapred <- function(data = NULL, model = NULL, label_known = NULL, amean = NULL, ascale = NULL, cmean = NULL, cscale = NULL) {
    raw_predict <- plspred(data, model, amean, ascale, cmean, cscale)
    raw_pred <- raw_predict$pred
    ns <- dim(raw_pred)[1]
    nYv <- dim(raw_pred)[2]
    label_predict <- matrix(0, ns, 1)    
    for (i in 1:ns) {
        label_predict[i] <- which(raw_pred[i,] == max(raw_pred[i,]))
    }
    results <- list(label_predict = label_predict, ccr = NULL, conf_mat = NULL)
    if (!is.null(label_known)) {
        ccr <- length(which(label_predict == label_known)) / length(label_known)
        conf_mat <- matrix(0, nYv, nYv)
        for (i in 1:nYv) {
            for (ii in 1:nYv) {
                conf_mat[i, ii] <- length(which(label_known == i & label_predict == ii)) / length(which(label_known == i))
            }
        }
        results$ccr <- ccr
        results$conf_mat <- conf_mat
    }
    return(results)
}



boots <- function(rep_idx = NULL) {
    library(pracma)
    rep_idx <- as.matrix(rep_idx)
    ns <- length(rep_idx)
    unique_rep <- unique(rep_idx)
    nr <- length(unique_rep)
    boots_idx <- unique(randi(nr, nr, 1))
    trn_idx <- which(!is.na(match(rep_idx, unique_rep[boots_idx])))
    tst_idx <- 1:ns
    tst_idx <- tst_idx[-trn_idx]
    results <- list(trn_idx = trn_idx, tst_idx = tst_idx)
    return(results)
}