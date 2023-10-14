#################################################################################################
#################################### GloWassRegfinal FUNCTION ###################################
#################################################################################################

#' @title Global Wasserstein Regression
#' @description  Global Frechet regression with respect to the Wasserstein distance.
#' @param xin An n by p matrix with input measurements of the predictors.
#' @param Qin An n by m matrix with values of quantile functions of which each row holds the 
#' quantile function values of an equally-spaced grid on [0, 1].
#' @param xout A k by p matrix with output measurements of the predictors.
#' @param lower A scalar with the lower bound of the support of the distribution. Default is
#'  \code{NULL}.
#' @param upper A scalar with the upper bound of the support of the distribution. Default is 
#' \code{NULL}.
#' @param Rsquared A logical variable indicating whether R squared would be returned. Default is 
#' FALSE.
#' @param Qgrid A numerical vector of length m holding the probability grid on [0, 1] at which
#'  the input quantile
#' functions take values. If \code{Rsquared} is TRUE, \code{Qgrid} is needed. Default is 
#' \code{seq(1,2*m,2)/2/m}.
#' @param weight survey weight train set

GloWassRegfinal = function(xin, Qin, xout, weight, lower = NULL, upper = NULL, Rsquared = TRUE,
                           Qgrid = NULL, plot = T, solver = T){
  require(fda.usc)
  if(is.matrix(Qin) == F){
    Qin = as.matrix(Qin)
  }
  
  if(is.data.frame(xin) == F){
    xin = as.data.frame(xin)
  }
  if(is.data.frame(xout) == F){
    xout = as.data.frame(xout)
  }
  if(nrow(xin) != nrow(Qin))
    stop("xin and Qin should have the same number of rows.")
  if(ncol(xin) != ncol(xout))
    stop("xin and xout should have the same number of columns.")
  if(Rsquared & is.null(Qgrid)){
    warning("Qgrid is missing and taking the default value.")
  }
  
  k = nrow(xout)
  n = nrow(xin) # Number of observations
  m = ncol(Qin) # Length of the quantile grid
  xbar = colMeans(xin) # Mean vector of the covariates
  Sigma = cov(xin) * (n - 1) / n # Covariance matrix
  invSigma = solve(Sigma)
  
  # If lower & upper are neither NULL
  A = cbind(diag(m), rep(0, m)) + cbind(rep(0, m), -diag(m)) # m + 1 columnas
  if(!is.null(upper) & !is.null(lower)){
    b0 = c(lower, rep(0, m - 1), -upper)
  }else if(!is.null(upper)){
    A = A[, -1]
    b0 = c(rep(0, m - 1), -upper)
  }else if(!is.null(lower)){
    A = A[, -ncol(A)]
    b0 = c(lower, rep(0, m - 1))
  }else{
    A = A[, -c(1, ncol(A))] # m - 1 columnas
    b0 = rep(0, m - 1)
  }
  
  Qout = sapply(1:k, function(j){
    s = 1 + (t(t(xin) - xbar)) %*% invSigma %*% t(xout[j, ] - xbar) # Frechet weights function
    s = as.vector(s)
    #' We standardised s*weight
    gx =  t(as.matrix(s*weight/sum(s*weight))) %*% Qin
    if(solver == T){
      res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
      return(sort(res$solution))
    }else{
      return(gx)
    }
  })
  
  Qout = t(Qout)
  # print(dim(Qout))
  
  if(!Rsquared) {
    return(Qout)
  }else{
    Qin.est = sapply(1:n, function(j){
      s = 1 + (t(t(xin) - xbar)) %*% invSigma %*% t(xin[j, ] - xbar)
      s = as.vector(s)
      gx =  t(as.matrix(s*weight/sum(s*weight))) %*% Qin
      if(solver == T){
        res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
        return(sort(res$solution))
      }else{
        return(gx)
      }
    })
    
    Qin.res =  t(Qin.est)
    
    frechetmedian = (weight/sum(weight)) %*% Qin.res
    if(plot){
      plot(fdata(frechetmedian))
    }
    
    errores2 = matrix(0, nrow = n, ncol = m)
    errores3 = matrix(0, nrow = n, ncol = m)
    
    for(i in 1:n){
      errores2[i, ] = (Qin[i, ] - frechetmedian) * (Qin[i, ] - frechetmedian)
    }
    
    for(i in 1:n){
      errores3[i, ] = (Qin[i, ] - Qin.res[i, ]) * (Qin[i, ] - Qin.res[i, ])
    }
    
    sumaerrores2 = apply(errores2, 1, mean)
    sumaerrores3 = apply(errores3, 1, mean)
    
    r2 = 1 - sum(weight/sum(weight) * sumaerrores3)/sum(weight/sum(weight) * sumaerrores2)
    return(list("Pred" = Qout, "r2" = r2, "Estim" = Qin.est))
  }
}

#################################################################################################
################################## concurrentconformal FUNCTION #################################
#################################################################################################

library(R.utils)
library(data.table)
library(tidyverse)
library(progress)

concurrente = function(formula, data, xpred = NULL, var_weights = rep(1, dim(data)[1]), 
                       silent = FALSE){
  ## Organize the input
  model_formula <- as.formula(formula)
  stopifnot(as.character(model_formula[[1]]) == '~' & length(model_formula) == 3)
  
  if(silent == FALSE){
    print("------------------------")
    print(model_formula)
  } 
  
  ##########################################################################################
  ## Step 1
  ##########################################################################################
  if(silent == FALSE) print("Step 1: Massive Univariate Model")
  
  covs <- all.vars(model_formula)[-1]
  X = data[, covs]
  if(length(covs) == 1){
    X <- as.data.frame(X)
    colnames(X) <- covs
  }
  var <- deparse(substitute(var_weights))
  tryCatch({
    if(is.numeric(var_weights) == T){
    }
  }, error = function(e){
    var_weights <<- data[[var]]
  })
  X$var_weights <- var_weights
  
  Y = data[,  all.vars(model_formula)[1]]
  
  
  L <- ncol(Y) ## number of observations on the functional domain
  argvals <- 1:L
  ## function "unimm" fit univariate mixed model at location l
  pos <- unlist(gregexpr('~', formula))
  formulas <- substr(formula, start = pos + 1, stop = nchar (formula))
  formulas <- paste0("Y ~", formulas)
  unimm <- function(l){
    data$Yl <- unclass(Y[, l])
    aux = data.frame("X" = X, "Y" = data$Yl)
    names(aux) <- sub("X.", "", names(aux))
    fit_uni <- lm(formula = formulas, data = aux, weights = var_weights)
    betaTilde <- fit_uni$coefficients
    residuals =  fit_uni$residuals
    rsquare <- summary(fit_uni)$adj.r.squared
    predict =  fit_uni$fitted.values
    if(is.null(xpred)){
      return(list("betatilde" = betaTilde, "residuals" = residuals, "r2" = rsquare,
                  "prediccion" = predict))
    }else{
      predict2 = predict(fit_uni, newdata = xpred, type = "response")
      # predict2 = as.matrix(predict2)
      # predict2 = as.data.frame(predict2)
      predict2 = as.numeric(unlist(predict2))
      names(predict2) = 1:length(predict2)
      
      return(list("betatilde" = betaTilde, "residuals" = residuals, "r2" = rsquare,
                  "prediccion" = predict, "prediccion2" = predict2))
    }
  }
  
  massmm = 1:length(argvals)
  massmm = as.list(massmm)
  for(i in 1:length(argvals)){
    massmm[[i]] = unimm(i)
    if(!(silent)) print(massmm[[i]])
  }
  # massmm <- lapply(argvals, unimm)
  
  ## obtain betaTilde
  require("dplyr")
  betaTilde <- t(lapply(massmm, '[[', 1) %>% bind_rows())
  r2 = t(lapply(massmm, '[', 3) %>% bind_rows())
  residuos =   t(lapply(massmm, '[[', 2) %>% bind_rows())
  predicciones =   t(lapply(massmm, '[[', 4) %>% bind_rows())
  colnames(betaTilde) <- argvals
  # print(betaTilde)
  
  if(is.null(xpred)){
    return(list("r2" = r2, "residuals" = residuos, "betatilde" = betaTilde,
                "prediccion" = predicciones))
  }else{
    predicciones2 = t(lapply(massmm, '[[', 5) %>% bind_rows())
    return(list("r2" = r2, "residuals" = residuos, "betatilde" = betaTilde,
                "prediccion" = predicciones, "prediccion2" = predicciones2))
  }
}

#################################################################################################

concurrentconformal = function(formula = "Y ~ X1 + X2 + X3 + X4 + X5", prob = 0.8, 
                               xpred = data[1, ], var_weights = rep(1, dim(data)[1]), 
                               data = data, silent = T, plot = T, tmin = NULL, tmax = NULL,
                               xlab = "Probability", ylab = "Estimated outcome",
                               idvar, id, solver = T, ...){
  model_formula <- as.formula(formula)
  stopifnot(as.character(model_formula[[1]]) == '~' & length(model_formula) == 3)
  var <- deparse(substitute(var_weights))
  tryCatch({
    if(is.numeric(var_weights) == T){
      if(sum(var_weights < 0) != 0) stop("Some weights are negative")
      # var_weights <- var_weights/sum(var_weights)
    }else if(var %in% names(data)){
      if(sum(data[[var]] < 0) != 0) stop("Some weights are negative")
      # data[[var]] <<- data[[var]]/sum(data[[var]])
      # var_weights <- data[[var]]/sum(data[[var]])
      var_weights <- data[[var]]
    }
  }, error = function(e){
    if(sum(data[[var]] < 0) != 0) stop("Some weights are negative")
    # data[[var]] <<- data[[var]]/sum(data[[var]])
    # var_weights <<- data[[var]]/sum(data[[var]])
    var_weights <<- data[[var]]
  })
  
  covs <- all.vars(model_formula)[-1]
  train <- data[, c(covs, all.vars(model_formula)[1])]
  rm(data)
  train$var_weights <- var_weights
  #' "They would be datasets "train", "cal1" and "cal2"
  xpred <- xpred[, covs]
  if(length(covs) == 1){
    xpred <- as.data.frame(xpred)
    colnames(xpred) <- covs
  }
  xpred[[all.vars(model_formula)[1]]] <- matrix(0, nrow = nrow(xpred), 
                                                ncol = ncol(train[[all.vars(model_formula)[1]]]))
  #' Splitting "train" into "train", "cal1" and "cal2"
  n1 <- ceiling(dim(train)[1]/3)
  n2 <- dim(train)[1] - n1
  n3 <- ceiling(n2/2)
  n2 <- n2 - n3
  # set.seed(seed)
  train$sample <- sample(c(rep(0, n1), rep(1, n2), rep(2, n3)), replace = F)
  cal1 <- train[train$sample == 1, ]
  cal1 <- cal1[, -dim(cal1)[2]]
  cal2 <- train[train$sample == 2, ]
  cal2 <- cal2[, -dim(cal2)[2]]
  var_weights1 <- var_weights[train$sample == 1]
  var_weights1 <- var_weights1/sum(var_weights1)
  var_weights2 <- var_weights[train$sample == 2]
  var_weights2 <- var_weights2/sum(var_weights2)
  var_weights <- var_weights[train$sample == 0]
  var_weights <- var_weights/sum(var_weights)
  train <- train[train$sample == 0, ]
  train <- train[, -dim(train)[2]]
  train$var_weights <- var_weights
  
  Y = train[,  all.vars(model_formula)[1]]
  # concurrentconf = concurrente(formula, xpred = cal1, data = train, var_weights = var_weights, 
  #                              silent = silent)
  require(fda.usc)
  g1 <- GloWassRegfinal(xin = train[, covs], Qin = Y, xout = cal1[, covs], weight = var_weights,
                        plot = F, solver = solver)
  media1 = g1$Pred
  Y1 = cal1[,  all.vars(model_formula)[1]]
  
  cal1[[all.vars(model_formula)[1]]] = abs(Y1 - media1)
  
  concurrentconf1 = concurrente(formula = formula, data = cal1, xpred = xpred,
                                var_weights = var_weights1, silent = silent)
  
  # sd = sqrt(concurrentconformal2$prediccion) ##### NA's (hai predici?ns negativas)
  sd = abs(concurrentconf1$prediccion2)
  
  
  g2 <- GloWassRegfinal(xin = train[, covs], Qin = Y, xout = cal2[, covs], weight = var_weights,
                        plot = F, solver = solver)
  media2 = g2$Pred
  Y2 = cal2[,  all.vars(model_formula)[1]]
  
  cal2[[all.vars(model_formula)[1]]] = abs(Y2 - media2)
  
  concurrentconf2 = concurrente(formula = formula, data = cal1, xpred = cal2, var_weights = var_weights1, 
                                silent = silent)
  
  sd2 = abs(concurrentconf2$prediccion2)
  est2 = (Y2 - media2)/sd2
  est2 = apply(est2, 1, max)
  prob <- 1 - (1 - prob)/2 # Para que deixe alpha/2 por riba e alpha/2 por baixo
  # cuantil = quantile(estandarizado2, probs = prob, na.rm = T)
  require("Hmisc")
  cuantil <- Hmisc::wtd.quantile(est2, var_weights2, prob = prob, na.rm = T, 
                                 normwt = TRUE)
  # bandasinqtrain = media - cuantil * sd
  # bandasdertrain = media + cuantil * sd
  
  g3 <- GloWassRegfinal(xin = train[, covs], Qin = Y, xout = xpred[, covs], weight = var_weights,
                        plot = F, solver = solver)
  
  media <- g3$Pred
  # sdtest = sqrt(concurrentconformal2$prediccion2)
  
  bandas1 = media - cuantil * sd
  bandas2 = media + cuantil * sd
  
  
  
  if(plot == T){
    if(is.null(tmin)){
      tmin <- 0
    }
    if(is.null(tmax)){
      tmax <- 1
    }
    t <- seq(tmin, tmax, length = dim(media)[2])
    if(missing(id)){
      id <- 1:dim(media)[1]
    }
    for(i in id){
      ylim = c(min(bandas1[i, ]), max(bandas2[i, ]))
      # ylim = c(0, 200)
      if(!missing(idvar)){
        main <- paste0("ID = ", xpred[[idvar]][i])
        plot(t, media[i, ], type = "l", ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
      }else{
        plot(t, media[i, ], type = "l", ylim = ylim, xlab = xlab, ylab = ylab)
      }
      polygon(c(t, rev(t)), c(bandas1[i, ], rev(bandas2[i, ])), border = NA, col = "lightgrey")
      lines(t, media[i, ], type = "l")
    }
  }
  
  res <- sqrt(mean((Y - t(g3$Estim))^2))
  return(list("Yest" = media, "upper" = bandas2, "lower" = bandas1, "r2" = g3$r2, 
              "quantile" = cuantil, "residuals" = res, "weights" = var_weights))
}
