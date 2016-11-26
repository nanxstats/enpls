#' Ensemble Partial Least Squares Regression
#'
#' Ensemble partial least squares regression.
#'
#' @param x Predictor matrix.
#' @param y Response vector.
#' @param maxcomp Maximum number of components included within each model.
#' If not specified, will use the maximum number possible (considering
#' cross-validation and special cases where n is smaller than p).
#' @param cvfolds Number of cross-validation folds used in each model
#' for automatic parameter selection, default is \code{5}.
#' @param reptimes Number of models to build with Monte-Carlo resampling
#' or bootstrapping.
#' @param method Resampling method. \code{"mc"} (Monte-Carlo resampling)
#' or \code{"boot"} (bootstrapping). Default is \code{"mc"}.
#' @param ratio Sampling ratio used when \code{method = "mc"}.
#' @param parallel Integer. Number of CPU cores to use.
#' Default is \code{1} (not parallelized).
#'
#' @return A list containing all partial least squares model objects.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enpls.fs}} for measuring feature importance
#' with ensemble partial least squares regressions.
#' See \code{\link{enpls.od}} for outlier detection with ensemble
#' partial least squares regressions.
#'
#' @export enpls.fit
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' fit = enpls.fit(x, y, reptimes = 50)
#' print(fit)
#' predict(fit, newx = x)

enpls.fit = function(x, y,
                     maxcomp  = NULL,
                     cvfolds  = 5L,
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
      plsdf = as.data.frame(cbind(xtmp, 'y' = ytmp))
      modellist[[i]] = suppressWarnings(enpls.fit.core(plsdf, maxcomp, cvfolds))
    }

  } else {

    registerDoParallel(parallel)
    modellist = foreach(i = 1L:reptimes) %dopar% {
      xtmp = x[samp.idx[[i]], ]
      ytmp = y[samp.idx[[i]]]
      plsdf = as.data.frame(cbind(xtmp, 'y' = ytmp))
      enpls.fit.core(plsdf, maxcomp, cvfolds)
    }

  }

  names(modellist) = paste0('pls_model_', 1L:length(modellist))

  class(modellist) = 'enpls.fit'
  return(modellist)

}

#' core function for enpls.fit
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set.
#' scale = TRUE
#'
#' @return the coefficients
#'
#' @keywords internal

enpls.fit.core = function(plsdf, maxcomp, cvfolds) {

  if (is.null(maxcomp)) {

    plsr.cvfit = plsr(y ~ .,
                      data       = plsdf,
                      scale      = TRUE,
                      method     = 'simpls',
                      validation = 'CV',
                      segments   = cvfolds)

  } else {

    plsr.cvfit = plsr(y ~ .,
                      data       = plsdf,
                      ncomp      = maxcomp,
                      scale      = TRUE,
                      method     = 'simpls',
                      validation = 'CV',
                      segments   = cvfolds)

  }

  # select best component number using adjusted CV
  cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])

  # remove plsr.cvfit object
  rm(plsr.cvfit)

  plsr.fit = plsr(y ~ .,
                  data       = plsdf,
                  ncomp      = cv.bestcomp,
                  scale      = TRUE,
                  method     = 'simpls',
                  validation = 'none')

  # minify plsr.fit object to reduce memory footprint
  plsr.fit[['model']] = NULL

  # save cv.bestcomp for predict.enpls
  enpls.core.fit = list('plsr.fit' = plsr.fit, 'cv.bestcomp' = cv.bestcomp)

  # remove plsr.fit object
  rm(plsr.fit)

  return(enpls.core.fit)

}
