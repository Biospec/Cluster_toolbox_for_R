# General use modelvalidation function 
#         modelvalidation <- function(train_fun = NULL, pred_fun = NULL, eval_type = 'R',
#                                    resample_method = 'boots', k= 1000, X = NULL, 
#                                    Y = NULL, rep_idx = NULL, model_parameters = NULL, 
#                                    cv_perm = FALSE, perm_test = FALSE)
#
# The name of training function is given in train_fun 
#  This function should follow a general format in which the first argument is the data matrix X and 
#  the second argument is the target vector/matrix Y. 
# The name of prediction function is given in pred_fun. 
#  This function should follow a general format in which first argument is the trained model generated 
#  by train_fun and the second argument is the testing data matrix
# 
# If eval_type is set to "R" (default) then a regression problem is assumed otherwise a classification problem 
#  is assumed. 
# For regression model, the function returns two performance statistics: Root mean square errors ($RMSEP) 
#  and coefficient of determination ($R2P)
# For classification model, the function returns correct classification rate ($CCR) and various classification related stats
#  ($stats)
#
# resample_method can be either "boots" for bootstrapping or "crossval" for crossvalidation
# If resample_method is "boots" then k is the number of iterations to be evaluated and set to 1000 by default.
# If resample_method is "crossval" then k is the number of folds cross-validation to be evaluated, if k = number of samples then
#  leave-one-out crossvalidation will be performed.
#
# X is the data matrix and Y is the target vector/matrix
# rep_idx is a vector showing which samples were the analytical replicates from the same sample. For example if five
#  samples were analysed and each sample had three repeated measurements then rep_idx =c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5).
#  If rep_idx is ignored then it is assumed that each sample had been analysed only once.
#
# cv_perm is to control whether perform a permutation before doing the cross-validation, if it is true the order of the samples
#  will be randomly shuffled so that the outcome would be less deterministic. Set to FALSE by default
#
# perm_test is a logical argument, if it is TRUE then a set of paired permutation tests will be also performed

library(caret)
library(permute)
library(pracma)

modelvalidation <- function(train_fun = NULL, pred_fun = NULL, eval_type = 'R',
                            resample_method = 'boots', k= 1000, X = NULL, 
                            Y = NULL, rep_idx = NULL, model_parameters = NULL, 
                            cv_perm = FALSE, perm_test = FALSE)
{ X <- as.matrix(X)
  ns <- dim(X)[1]
# Only handle numerical Y at the moment
#  if (is.numeric(Y)){
    Y <- as.matrix(Y)
    nvY <- dim(Y)[2]
    if (nvY > 1 & eval_type != 'R'){
      warning("Classification target Y appeared to be a matrix, convert it to vector by taking maximum on each row",
              immediate = TRUE)
      Ynew = matrix(0, ns, 1)
      for (i in 1:ns) Ynew[i,] <- which(Y[i,] == max(Y[i,]))
      Y <- Ynew
      nvY <- 1
    }
#  }
  # For regression problems, each unique concentration or concentrations combination would be considered as a rep.
  if (eval_type == "R" & nvY == 1)
    rep_idx <- Y
  if (eval_type == "R" & nvY > 1){
    unique_Y <- unique(Y)
    nuY <- dim(unique_Y)[1]
    rep_idx <- matrix(0,nuY,1)
    for (i in 1:nuY){
      nuY_idx <- which(data.frame(t(Y)) %in% data.frame(t[unique_Y[i,]]))
      rep_idx[nuY_idx] <- i
    }
  }
  
  if (is.null(rep_idx)) 
    rep_idx = 1:ns
  
  if (resample_method == 'boots'){
    if (eval_type == 'R') {
      rmsep_boots <- matrix(0, k, nvY)
      r2p_boots <- matrix(0, k, nvY)
    } else {
      ccr <- matrix(0, k, 1)
      confmat <- list()
    }
    
    boots_idx <- boots(rep_idx)
    for (i in 1:k){
      trn_idx <- boots_idx$trn_idx[[i]]
      tst_idx <- boots_idx$tst_idx[[i]]
      opt_model <- modeltune(train_fun, pred_fun, eval_type, 
                                  k, X[trn_idx,], Y[trn_idx,], rep_idx[trn_idx,], 
                                  model_parameters = model_parameters)
      opt_parameters <- opt_model$opt_parameters
      model <- do.call(train_fun, c(X[trn_idx,], Y[trn_idx,], opt_parameters))
      pred <- pred_fun(model, X[tst_idx,])
      
      if (eval_type == 'R') {
        y_known <- Y[boots_idx$tst_idx,]
        y_predict <- pred
        msep_boots[i,] <- sum((y_pred - y_known)^2)
        r2p_boots[i,] <- rmsep^2 / sum((y_known - colMeans(y_known))^2) 
        rmsep_boots <- sqrt(msep_boots)
      } else{
        y_known <- Y[boots_idx$tst_idx,]
        y_predict <- pred
        ccr[i] <- length(which(y_predict == y_known))/length(y_known)
        confmat <- c(confmat, confusionMatrix(y_predict, y_known))
      }
    }
    
    if (eval_type == 'R') 
      results <- c(RMSEP = rmsep_boots, R2P = r2p_boots)
    else
      results <- c(CCR = ccr, ConfMat = confmat)
  }
  else if (resample_method == 'crossval'){
    cv_idx <- crossval(rep_idx, k, cv_perm)
    y_pred_cv <- matrix(0, ns, nvY)
    y_known_cv <- matrix(0, ns, nvY)
    no_loops <- length(cv_idx$trn_idx)
    for (i in 1:no_loops){
      trn_idx <- cv_idx$trn_idx[[i]]
      tst_idx <- cv_idx$tst_idx[[i]]
      model <- train_fun(X[trn_idx,], Y[trn_idx,],...)
      pred <- pred_fun(model, X[tst_idx,])
      y_known_cv[tst_idx,] <- Y[tst_idx,]
      y_pred_cv[tst_idx,] <- pred
    }
    
    if (eval_type == 'R'){
      msep_cv <- sum((y_pred_cv - y_known_cv)^2) / ns
      r2cv <- msep_cv / sum((y_knwon - colMeans(y_known_cv))^2)
      rmsep_cv <- sqrt(msep_cv)
    } else {
      ccr <- length(which(y_pred_cv == y_known_cv))/length(y_known_cv)
      confmat <- confusionMatrix(y_pred_cv, y_known_cv)
    }
    
    if (eval_type == 'R')
      results <- c(RMSEP = rmsep_cv, R2P = r2cv)
    else
      results <- c(CCR = ccr, ConfMat = confmat)
  } else 
    stop("The resample_method can only be either boots or crossval")
  return(results)
}

# Model selection function via cross-validation, arguments are the same as modelvalidation

modeltune <- function(train_fun = NULL, pred_fun = NULL, eval_type = "R", 
                      k = 7, X = NULL, Y = NULL, rep_idx = NULL, 
                      cv_perm = FALSE, perm_test = FALSE, model_parameters = NULL){
  X <- as.matrix(X)
  ns <- dim(X)[1]
# Only handle numerical Y at the moment  
#  if (is.numeric(Y) | is.data.frame(Y)){
    Y <- as.matrix(Y)
    nvY <- dim(Y)[2]
    if (nvY > 1 & eval_type != 'R'){
      warning("Classification target Y appeared to be a matrix, convert it to vector by taking maximum on each row",
              immediate = TRUE)
      Ynew = matrix(0, ns, 1)
      for (i in 1:ns) Ynew[i,] <- which(Y[i,] == max(Y[i,]))
      Y <- Ynew
      nvY <- 1
    }
#  }
  # For regression problems, each unique concentration or concentrations combination would be considered as a rep.
  if (eval_type == "R" & nvY == 1)
    rep_idx <- Y
  if (eval_type == "R" & nvY > 1){
    unique_Y <- unique(Y)
    nuY <- dim(unique_Y)[1]
    rep_idx <- matrix(0,nuY,1)
    for (i in 1:nuY){
      nuY_idx <- which(match(data.frame(t(unique_Y[i,])), data.frame(t[Y])) == 1)
      rep_idx[nuY_idx] <- i
    }
  }
    
  no_parameters <- length(model_parameters)
  if (no_parameters > 1) {
    parameter_grid <- expand.grid(model_parameters)
    no_trials <- dim(parameter_grid)[1]
  }
  else {
    parameter_grid <- model_parameters
    no_trials <- length(parameter_grid[[1]])
  }
  
  if (is.null(rep_idx)) 
    rep_idx <- 1:ns
  
  if (eval_type == 'R') {
    r2cv <- matrix(0, no_trials, nvY)
    rmsep_cv <- matrix(0, no_trials, nvY)
  } else {
    ccr <- matrix(0, k, 1)
    clsstats <- list()
  }
  
  cv_idx <- crossval(rep_idx, k, cv_perm)
  no_loops <- length(cv_idx$trn_idx)
  
  for (i in 1:no_trials){
    y_pred_cv <- matrix(0, ns, nvY)
    y_known_cv <- matrix(0, ns, nvY)
    no_loops <- length(cv_idx$trn_idx)
    for (ii in 1:no_loops){
      trn_idx <- unlist(cv_idx$trn_idx[[ii]])
      tst_idx <- unlist(cv_idx$tst_idx[[ii]])
      if (is.data.frame(parameter_grid)){
        model <- do.call(train_fun, c(list(X[trn_idx, ], Y[trn_idx, ]), as.list(parameter_grid[i,])))
      }
      else{
        par_name <- names(parameter_grid)

        par_val <- parameter_grid[[1]]
        par_list <- eval(parse(text = paste('c(list(X[trn_idx, ], Y[trn_idx, ]), ', 'list(', par_name, ' = ', par_val[i], '))', sep = "")))
        model <- do.call(train_fun, par_list)
      }
      pred <- pred_fun(model, X[tst_idx, ])
      y_known_cv[tst_idx, ] <- Y[tst_idx, ]
      y_pred_cv[tst_idx, ] <- pred
    }
    if (eval_type == 'R'){
      msep_cv <- sum((y_pred_cv - y_known_cv)^2) / ns
      if (nvY > 1) {
        r2cv[i,] <- 1 - msep_cv / (sum((y_knwon_cv - colMeans(y_known_cv))^2) / ns)
        rmsep_cv[i,] <- sqrt(msep_cv)
      } else {
        r2cv[i] <- 1 - msep_cv / (sum((y_known_cv - mean(y_known_cv))^2) / ns)
        rmsep_cv[i] <- sqrt(msep_cv)
      }
    } else {
      ccr[i] <- length(which(y_pred_cv == y_known_cv))/length(y_known_cv)
      clsstats <- c(clsstats, confusionMatrix(y_pred_cv, y_known_cv))
    }
  }
  if (eval_type == 'R'){
    if (nvY == 1) {
      opt_idx <- which(r2cv==max(r2cv))
      if (length(opt_idx) > 1) 
        opt_idx <- opt_idx[1]
      opt_results <- c(r2cv = r2cv[opt_idx], rmsecv = rmsep_cv[opt_idx])
    } else {
      r2cv_mean <- rowMeans(r2cv)
      opt_idx <- which(r2cv_mean == max(r2cv))
      if (length(opt_idx) > 1)
        opt_idx <- opt_idx[1]
      opt_results <- c(r2cv = r2cv[opt_idx,], rmsecv = rmsep_cv[opt_idx,])
    }
  } else {
    opt_idx <- which(ccr == max(ccr))
    if (length(opt_idx) > 1)
      opt_idx <- opt_idx[1]
    opt_results <- list(CCR = ccr[opt_idx], stats = clsstats[opt_idx])
  }
  if (no_parameters > 1)
    results <- list(opt_parameters = as.list(parameter_grid[opt_idx, ]), opt_results = opt_results)
  else {
    opt_parameters <- eval(parse(text = paste('list(', par_name, ' = ', par_val[opt_idx], ')', sep = '')))
    results <- list(opt_parameters = opt_parameters, opt_results = opt_results)
  }
  
  return(results)
}
  
# Generate subsets of samples via bootstrapping
boots <- function(rep_idx = NULL, k = 1000) {
  rep_idx <- as.matrix(rep_idx)
  ns <- length(rep_idx)
  unique_rep <- unique(rep_idx)
  nr <- length(unique_rep)
  trn_idx_boots <- list()
  tst_idx_boots <- list()
  for (i in 1:k){
    boots_idx <- randi(nr, nr, 1)
    trn_idx <- list()
    sample_idx <- 1:ns
    for (ii in 1:nr)
      trn_idx[[ii]] <- sample_idx[which(match(rep_idx, unique_rep[boots_idx[ii]])==1)]
    trn_idx <- sort(unlist(trn_idx))
    tst_idx <- sample_idx[-trn_idx]
    trn_idx_boots[[i]] <- trn_idx
    tst_idx_boots[[i]] <- tst_idx
  }
  results <- list(trn_idx = trn_idx_boots, tst_idx = tst_idx_boots)
  return(results)
}

# Generate subsets of samples via k-cross-validation
crossval <- function(rep_idx = NULL, k = NULL, perm = FALSE) {
  rep_idx <- as.matrix(rep_idx)
  ns <- length(rep_idx)
  unique_rep <- unique(rep_idx)
  nr <- length(unique_rep)
  if (is.null(k) | k > nr) k <- nr
  if (perm) {
    perm_idx <- shuffle(1:nr)
    unique_rep <- unique_rep[perm_idx]
  }
  step_size <- round(nr / k)
  trn_idx_cv <-list()
  tst_idx_cv <-list()
  count <- 1
  
  for (i in seq(1, nr, step_size)){
    rep_cv <- unique_rep[i:min(nr, i+step_size-1)]
    tst_idx <- which(!is.na(match(rep_idx, rep_cv)))
    trn_idx <- 1:ns
    trn_idx <- trn_idx[-tst_idx]
    trn_idx_cv[[count]] <- as.list(trn_idx)
    tst_idx_cv[[count]] <- as.list(tst_idx)
    count <- count + 1
  }

  results <- list(trn_idx = trn_idx_cv, tst_idx = tst_idx_cv)
  return(results)
}

# Rep based permutation
reps_perm <- function(reps = NULL, Y = NULL) {
# Only handle numerical Y at the moment
#  if (is.numeric(Y))  {
      Y <- as.matrix(Y)
      Y_perm = matrix(0, dim(Y))
#  }
  
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
    Y_new <- Y[idx:(idx + no_reps_found - 1), ]
    if (dim(unique(Y_new))[1] > 1) Y_new[1:no_reps_found,] <- Y_new[1,] # ensure replicates after permutation still has the same label
    Y_perm[idx:(idx + no_reps_found - 1), ] <- Y_new
    idx <- idx + no_reps_found
  }
  
  return(c(perm_idx = perm_idx, Y_perm = Y_perm))
}
