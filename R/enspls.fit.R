#' Ensemble Sparse Partial Least Squares Regression
#'
#' Ensemble sparse partial least squares regression.
#'
#' @param x Predictor matrix.
#' @param y Response vector.
#' @param maxcomp Maximum number of components included within each model.
#' If not specified, will use \code{5} by default.
#' @param cvfolds Number of cross-validation folds used in each model
#' for automatic parameter selection, default is \code{5}.
#' @param alpha Parameter (grid) controlling sparsity of the model.
#' If not specified, default is \code{seq(0.2, 0.8, 0.2)}.
#' @param reptimes Number of models to build with Monte-Carlo resampling
#' or bootstrapping.
#' @param method Resampling method. \code{"mc"} (Monte-Carlo resampling)
#' or \code{"boot"} (bootstrapping). Default is \code{"mc"}.
#' @param ratio Sampling ratio used when \code{method = "mc"}.
#' @param parallel Integer. Number of CPU cores to use.
#' Default is \code{1} (not parallelized).
#'
#' @return A list containing all sparse partial least squares model objects.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{enspls.fs}} for measuring feature importance
#' with ensemble sparse partial least squares regressions.
#' See \code{\link{enspls.od}} for outlier detection with ensemble
#' sparse partial least squares regressions.
#'
#' @export enspls.fit
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @examples
#' data("logd1k")
#' x = logd1k$x
#' y = logd1k$y
#'
#' set.seed(42)
#' fit = enspls.fit(x, y, reptimes = 5, maxcomp = 3,
#'                  alpha = c(0.3, 0.6, 0.9))
#' print(fit)
#' predict(fit, newx = x)

enspls.fit = function(x, y,
                      maxcomp  = 5L,
                      cvfolds  = 5L,
                      alpha    = seq(0.2, 0.8, 0.2),
                      reptimes = 500L,
                      method   = c('mc', 'boot'),
                      ratio    = 0.8,
                      parallel = 1L) {

  if (missing(x) | missing(y)) stop('Please specify both x and y')

  method = match.arg(method)

  x.row = nrow(x)
  samp.idx = vector('list', reptimes)

  if (method == 'mc') {
    for (i in 1L:reptimes) samp.idx[[i]] = sample(1L:x.row, round(x.row * ratio))
  }

  if (method == 'boot') {
    for (i in 1L:reptimes) samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
  }

  if (parallel < 1.5) {

    modellist = vector('list', reptimes)
    for (i in 1L:reptimes) {
      xtmp = x[samp.idx[[i]], ]
      ytmp = y[samp.idx[[i]]]
      modellist[[i]] = enspls.fit.core(xtmp, ytmp, maxcomp, cvfolds, alpha)
    }

  } else {

    registerDoParallel(parallel)
    modellist = foreach(i = 1L:reptimes) %dopar% {
      xtmp = x[samp.idx[[i]], ]
      ytmp = y[samp.idx[[i]]]
      enspls.fit.core(xtmp, ytmp, maxcomp, cvfolds, alpha)
    }

  }

  names(modellist) = paste0('spls_model_', 1L:length(modellist))

  class(modellist) = 'enspls.fit'
  return(modellist)

}

#' core function for enspls.fit
#'
#' select the best ncomp and alpha, then use them to fit
#' the complete training set.
#' scale = TRUE
#'
#' @importFrom spls cv.spls spls
#' @importFrom utils capture.output
#'
#' @return the coefficients
#'
#' @keywords internal

enspls.fit.core = function(xtmp, ytmp, maxcomp, cvfolds, alpha) {

  invisible(
    capture.output(
      spls.cvfit <- cv.spls(xtmp,
                            ytmp,
                            fold    = cvfolds,
                            K       = maxcomp,
                            eta     = alpha,
                            scale.x = TRUE,
                            scale.y = FALSE,
                            plot.it = FALSE)))

  # select best component number and alpha using adjusted CV
  cv.bestcomp  = spls.cvfit$'K.opt'
  cv.bestalpha = spls.cvfit$'eta.opt'

  # clean up spls.cvfit object
  rm(spls.cvfit)

  spls.fit = spls(xtmp,
                  ytmp,
                  K       = cv.bestcomp,
                  eta     = cv.bestalpha,
                  scale.x = TRUE,
                  scale.y = FALSE)

  # save cv.bestcomp and cv.bestalpha for predict.enspls
  enspls.core.fit = list('spls.fit' = spls.fit,
                         'cv.bestcomp' = cv.bestcomp,
                         'cv.bestalpha' = cv.bestalpha)

  # clean up spls.fit object
  rm(spls.fit)

  return(enspls.core.fit)

}
