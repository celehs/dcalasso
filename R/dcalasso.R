#' Divide-and-conquer method for the fitting of adaptive lasso model with big data
#'
#' \code{dcalasso} fits adaptive lasso for big datasets using multiple linearization methods,
#' including one-step estimation and least square approximation. This function is able to
#' fit the adaptive lasso model either when the dataset is being loaded as a whole into \code{data} or when
#' the datasets are splitted a priori and saved into multiple \code{rds} files.
#' The algorithm uses a divide-and-conquer one-step estimator as the initial estimator
#' and uses a least square approximation to the partial likelihood, which
#' reduces the computation cost. The algorithm currently supports adaptive lasso with
#' Cox proportional hazards model with or without
#' time-dependent covariates. Ties in survival data analysis are handled by Efron's method.
#' The first half of the routine computes an initial estimator (n^{1/2} consistent estimator). It first obtains a warm-start by
#' fitting coxph to the first subset (first random split of data or first data file indicated by data.rds) and then uses one-step
#' estimation with iter.os rounds to update the warm-start. The one-step estimation loops through each subset and gathering scores
#' and information matrices. The second half of the routine then shrinks the initial estimator using a least square approximation-based adaptive lasso step.
#'
#' @param formula a formula specifying the model. For Cox model, the outcome should be specified as the Surv(start, stop, status) or Surv(start, status) object in the survival package.
#' @param family For Cox model, family should be cox.ph(), or "cox.ph".
#' @param data data frame containing all variables.
#' @param data.rds when the dataset is too big to load as a whole into the RAM, one can specify \code{data.rds}
#' which are the full paths of all randomly splitted subsets of the full data, saved into multiple \code{.rds} format.
#' @param weights a prior weights on each observation
#' @param subset an expression indicating subset of rows of data used in model fitting
#' @param na.action how to handle NA
#' @param offset an offset term with a fixed coefficient of one
#' @param lambda tuning parameter for the adaptive lasso penalty. penalty = lambda * sum_j |beta_j|/|beta_j initial|^gamma
#' @param gamma exponent of the adaptive penalty. penalty = lambda * sum_j |beta_j|/|beta_j initial|^gamma
#' @param K number of division of the full dataset. It will be overwritten to \code{length(data.rds)} if data.rds is given.
#' @param iter.os number of iterations for one-step updates
#' @param ncores number of cores to use. The iterations will be paralleled using \code{foreach} if ncores>1.
#' @return
#' \item{coefficients.pen}{adaptive lasso shrinkage estimation}
#' \item{coefficients.unpen}{initial unregularized estimator}
#' \item{cov.unpen}{variance-covariance matrix of unpenalized model}
#' \item{cov.pen}{variance-covariance matrix of penalized model}
#' \item{BIC}{sequence of BIC evaluation at each lambda}
#' \item{n.pen}{number use to penalize the degrees of freedom in BIC. }
#' \item{n}{number of used rows of the data}
#' \item{idx.opt}{index for the optimal BIC}
#' \item{BIC.opt}{minimal BIC}
#' \item{family}{family object of the model}
#' \item{lamba.opt}{optimal lambda to minimize BIC}
#' \item{df}{degrees of freedom at each lambda}
#' \item{p}{number of covariates}
#' \item{iter}{number of one-step iterations}
#' \item{Terms}{term object of the model}
#' @author Yan Wang \email{yaw719@mail.harvard.edu}, Tianxi Cai \email{tcai@hsph.harvard.edu}, Chuan Hong <Chuan_Hong@hms.harvard.edu>
#' @references Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "A Fast Divide-and-Conquer Sparse Cox Regression." arXiv preprint arXiv:1804.00735 (2018).
#' @export
dcalasso <- function(formula, family=cox.ph(), data = NULL, data.rds = NULL, weights, subset, na.action,
                     offset, lambda = 10^seq(-10,3,0.01), gamma = 1, K = 20, iter.os = 2, ncores = 1){

  if (is.null(data) & is.null(data.rds)){
    stop('Either data or data.rds needs to be specified.')
  }
  K = as.integer(K); iter.os = as.integer(iter.os)
  if (iter.os < 0){
    stop('iter.os<0 is not allowed.')
  }
  if (iter.os ==0 & K!=1){
    stop('iter.os=0 and K>1 is not allowed.')
  }
  if (is.character(family)){
    family = eval(parse(text = family))
  }
  if (is.function(family)){
    family = family()
  }
  if(is.null(family$family)){
    stop('family not recognized')
  }
  if (family$family != 'Cox PH')
      stop("The current version supports only Cox model.")
  family <- fix.family.link(family)
  family <- fix.family.var(family)
  family <- fix.family.ls(family)

  mf = match.call()
  mf$family <- mf$K <- mf$iter.os <- mf$gamma <- mf$lambda <- mf$ncores <- mf$data.rds <- NULL

  special <- c("strata")
  mf$formula <- if (missing(data)) terms(formula,special)
  else terms(formula, special, data = data)
  mf[[1]] = quote(stats::model.frame)

  if (!is.null(data)){
    ############# Dataset loaded as a whole ######################
    mf = eval(mf, parent.frame())
    rm(data)
    Terms = terms(mf)

    tmp = dataprocess(mf, Terms, family)
    Y = tmp$Y; X = tmp$X; strats = tmp$strats; weights = tmp$weights; offset = tmp$offset
    rm(mf); rm(tmp)

    # Warm-start
    options(warn = -1)
    set.seed(32767)
    dc.idx = split(1:dim(X)[1], factor(1:K))
    options(warn = 1)

    bini = warmfit(Y, X, strats, weights, offset, family, dc.idx[[1]])
    if (K == 1){
      tmp = scoreinfofun (Y, X, strats, weights, offset, family, bini, dc.idx[[1]])
      score = tmp$score; info = tmp$info
    }else{
      # Followed by one-step updates
      if (ncores == 1){
        for (iter in 1:iter.os){
          # Gather score and information
          score = rep(0, length(bini)); info = matrix(0, length(bini), length(bini))
          for (i in 1:length(dc.idx)){
            tmp = scoreinfofun (Y, X, strats, weights, offset, family, bini, dc.idx[[i]])
            score = score + tmp$score
            info = info + tmp$info
          }
          bini = as.vector(bini + solve(info) %*% score)
        }
      }else{
        cl = parallel::makeCluster(ncores)
        registerDoParallel(cl)
        for (iter in 1:iter.os){
          # Gather score and information
          tmp.combine = foreach (i = 1:length(dc.idx), .combine = rbind) %dopar%{
            tmp = scoreinfofun (Y, X, strats, weights, offset, family, bini, dc.idx[[i]])
            return(list(score = tmp$score, info = tmp$info))
          }
          score = rowSums(sapply(1:length(dc.idx), function(kk) {
            tmp.combine[[kk,1]]
          }))
          info = Reduce('+',lapply(1:length(dc.idx), function(kk) {
            tmp.combine[[kk,2]]
          }))
          bini = as.vector(bini + solve(info) %*% score)
        }
        parallel::stopCluster(cl)
      }
    }
    if (family$family == 'Cox PH'){
      if (ncol(Y) == 2){
        n.pen = sum(Y[,2])
      }else{
        n.pen = sum(Y[,3])
      }
    }else{
      stop("The current version supports only Cox model.")
    }
    n = dim(X)[1]
  }else{
    ########### Datasets saved in separate files #######################
    mf0 = mf
    K = length(data.rds)
    if (ncores == 1){
      for (iter in 1:iter.os){
        for (i in 1:K){
          mf0$data = as.data.frame(readRDS(data.rds[i]))
          mf = eval(mf0, parent.frame())
          Terms = terms(mf)

          tmp = dataprocess(mf, Terms, family)
          Y = tmp$Y; X = tmp$X; strats = tmp$strats; weights = tmp$weights; offset = tmp$offset
          rm(mf); rm(tmp)

          if (i == 1){
            bini = warmfit(Y, X, strats, weights, offset, family, 1:length(weights))
            score = rep(0, length(bini)); info = matrix(0, length(bini), length(bini)); n.pen = n = 0
          }
          tmp = scoreinfofun (Y, X, strats, weights, offset, family, bini, 1:length(weights))
          score = score + tmp$score
          info = info + tmp$info
          if (family$family == 'Cox PH'){
            if (ncol(Y) == 2){
              n.pen = n.pen + sum(Y[,2])
            }else{
              n.pen = n.pen + sum(Y[,3])
            }
          }else{
            stop("The current version supports only Cox model.")
          }
          n = n + dim(X)[1]
        }
        bini = as.vector(bini + solve(info) %*% score)
      }
    }else{
      mf0$data = as.data.frame(readRDS(data.rds[1]))
      mf = eval(mf0, parent.frame())
      Terms = terms(mf)

      tmp = dataprocess(mf, Terms, family)
      Y = tmp$Y; X = tmp$X; strats = tmp$strats; weights = tmp$weights; offset = tmp$offset
      rm(mf); rm(tmp)

      bini = warmfit(Y, X, strats, weights, offset, family, 1:length(weights))
      cl = parallel::makeCluster(ncores)
      registerDoParallel(cl)
      for (iter in 1:iter.os){
        tmp.combine = foreach (i = 1:K, .combine = rbind) %dopar%{
          mf0$data = as.data.frame(readRDS(data.rds[i]))
          mf = eval(mf0, parent.frame())
          Terms = terms(mf)

          tmp = dataprocess(mf, Terms, family)
          Y = tmp$Y; X = tmp$X; strats = tmp$strats; weights = tmp$weights; offset = tmp$offset
          rm(mf); rm(tmp)

          tmp = scoreinfofun (Y, X, strats, weights, offset, family, bini, 1:length(weights))
          if (family$family == 'Cox PH'){
            if (ncol(Y) == 2){
              n.pen = sum(Y[,2])
            }else{
              n.pen = sum(Y[,3])
            }
          }else{
            stop("The current version supports only Cox model.")
          }
          n = dim(X)[1]
          return(list(score = tmp$score, info = tmp$info, n.pen = n.pen, n = n))
        }
        score = rowSums(sapply(1:K, function(kk) {
          tmp.combine[[kk,1]]
        }))
        info = Reduce('+',lapply(1:K, function(kk) {
          tmp.combine[[kk,2]]
        }))
        n.pen = sum(sapply(1:K, function(kk) {
          tmp.combine[[kk,3]]
        }))
        n = sum(sapply(1:K, function(kk) {
          tmp.combine[[kk,4]]
        }))
        bini = as.vector(bini + solve(info) %*% score)
      }
      parallel::stopCluster(cl)
    }
  }

  ######## Penalized regression
  infohalf = svd(info/n); infohalf = infohalf$u%*%diag(sqrt(infohalf$d))%*%t(infohalf$v);
  tmpfit = glmnet(y = infohalf %*% bini,
                  x = infohalf,
                  family='gaussian',
                  penalty.factor = 1/(abs(bini)^gamma),
                  alpha=1,
                  lambda = lambda,intercept=F)
  LL = apply((c(infohalf %*% bini) - predict(tmpfit,infohalf,type = 'response'))^2, 2, sum)*n
  # BIC
  BIC.lam = LL+log(n.pen)*tmpfit$df
  m.opt = which.min(BIC.lam); bhat.BIC = tmpfit$beta[,m.opt]; lamhat.BIC = tmpfit$lambda[m.opt]

  names(bini) = names(bhat.BIC) = colnames(info) = row.names(info) = colnames(X)
  cov.unpen = solve(info)
  cov.pen   = cov.unpen; cov.pen[bhat.BIC==0, ] = 0; cov.pen[,bhat.BIC==0] = 0;
  object = list(n.pen = n.pen, n = n, K = K,
                coefficients.unpen = bini, coefficients.pen = bhat.BIC,
                cov.unpen = cov.unpen, cov.pen = cov.pen, lambda = lambda,
                BIC = BIC.lam, idx.opt = m.opt, BIC.opt = BIC.lam[m.opt], family = family,
                lambda.opt = lamhat.BIC, df = tmpfit$df, p = length(bini), iter = iter.os,
                Terms = Terms)
  class(object) = 'dcalasso'
  return(object)
}
