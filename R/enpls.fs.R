#' Ensemble Partial Least Squares for Measuring Feature Importance
#'
#' Measuring feature importance with ensemble partial least squares.
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
#' @return A list containing two components:
#' \itemize{
#' \item \code{variable.importance} - a vector of variable importance
#' \item \code{coefficient.matrix} - original coefficient matrix
#' }
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enpls.od}} for outlier detection with
#' ensemble partial least squares regressions.
#' See \code{\link{enpls.fit}} for fitting ensemble partial least
#' squares regression models.
#'
#' @export enpls.fs
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
#' fs = enpls.fs(x, y, reptimes = 50)
#' print(fs)
#' plot(fs)

enpls.fs = function(x, y,
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

    coeflist = vector('list', reptimes)
    for (i in 1L:reptimes) {
      xtmp = x[samp.idx[[i]], ]
      xtmp = scale(xtmp, center = TRUE, scale = TRUE)
      ytmp = y[samp.idx[[i]]]
      plsdf = as.data.frame(cbind(xtmp, 'y' = ytmp))
      coeflist[[i]] = suppressWarnings(enpls.fs.core(plsdf, maxcomp, cvfolds))
    }

  } else {

    registerDoParallel(parallel)
    coeflist = foreach(i = 1L:reptimes) %dopar% {
      xtmp = x[samp.idx[[i]], ]
      xtmp = scale(xtmp, center = TRUE, scale = TRUE)
      ytmp = y[samp.idx[[i]]]
      plsdf = as.data.frame(cbind(xtmp, 'y' = ytmp))
      enpls.fs.core(plsdf, maxcomp, cvfolds)
    }

  }

  coefmat = do.call(rbind, coeflist)

  varimp = abs(colMeans(coefmat))/apply(coefmat, 2L, sd)

  object = list('variable.importance' = varimp,
                'coefficient.matrix'  = coefmat)
  class(object) = 'enpls.fs'
  return(object)

}

#' core function for enpls.fs
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set.
#' scale = FALSE
#'
#' @return fitted coefficients
#'
#' @keywords internal

enpls.fs.core = function(plsdf, maxcomp, cvfolds) {

  if (is.null(maxcomp)) {

    plsr.cvfit = plsr(y ~ .,
                      data       = plsdf,
                      scale      = FALSE,
                      method     = 'simpls',
                      validation = 'CV',
                      segments   = cvfolds)

  } else {

    plsr.cvfit = plsr(y ~ .,
                      data       = plsdf,
                      ncomp      = maxcomp,
                      scale      = FALSE,
                      method     = 'simpls',
                      validation = 'CV',
                      segments   = cvfolds)

  }

  # select best component number using adjusted CV
  cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])

  plsr.fit = plsr(y ~ .,
                  data       = plsdf,
                  ncomp      = cv.bestcomp,
                  scale      = FALSE,
                  method     = 'simpls',
                  validation = 'none')

  plsr.coef = drop(coef(plsr.fit))

  return(plsr.coef)

}
