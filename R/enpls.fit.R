#' Ensemble Partial Least Squares Regression
#'
#' Ensemble partial least squares regression.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param maxcomp Maximum number of components included within the models,
#' if not specified, default is the variable (column) numbers in x.
#' @param MCtimes times of Monte-Carlo
#' @param method \code{"mc"} or \code{"bootstrap"}. Default is \code{"mc"}.
#' @param ratio sample ratio used when \code{method = "mc"}
#' @param parallel Integer. Number of CPU cores to use.
#' Default is \code{1} (not parallelized).
#'
#' @return A list containing all partial least squares model objects.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble
#' partial least squares regression.
#' See \code{\link{enpls.od}} for outlier detection with ensemble
#' partial least squares regression.
#'
#' @export enpls.fit
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @references
#' Dongsheng Cao, Yizeng Liang, Qingsong Xu, Yifeng Yun, and Hongdong Li.
#' "Toward better QSAR/QSPR modeling: simultaneous outlier detection and
#' variable selection using distribution of model features."
#' \emph{Journal of computer-aided molecular design} 25, no. 1 (2011): 67--80.
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' fit = enpls.fit(x, y, MCtimes = 100)
#' print(fit)
#' predict(fit, newx = x)

enpls.fit = function(x, y,
                     maxcomp = NULL,
                     MCtimes = 500L,
                     method = c('mc', 'bootstrap'), ratio = 0.8,
                     parallel = 1L) {

  if (missing(x) | missing(y)) stop('Please specify both x and y')

  if (is.null(maxcomp)) maxcomp = ncol(x)

  method = match.arg(method)

  x.row = nrow(x)
  samp.idx = vector('list', MCtimes)

  if (method == 'mc') {
    for (i in 1L:MCtimes) samp.idx[[i]] = sample(1L:x.row, round(x.row * ratio))
  }

  if (method == 'bootstrap') {
    for (i in 1L:MCtimes) samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
  }

  if (parallel < 1.5) {

    modellist = vector('list', MCtimes)
    for (i in 1L:MCtimes) {
      xtmp = x[samp.idx[[i]], ]
      ytmp = y[samp.idx[[i]]]
      plsdf = as.data.frame(cbind(xtmp, 'y' = ytmp))
      modellist[[i]] = suppressWarnings(enpls.fit.core(plsdf, maxcomp))
    }

  } else {

    registerDoParallel(parallel)
    modellist = foreach(i = 1L:MCtimes) %dopar% {
      xtmp = x[samp.idx[[i]], ]
      ytmp = y[samp.idx[[i]]]
      plsdf = as.data.frame(cbind(xtmp, 'y' = ytmp))
      enpls.fit.core(plsdf, maxcomp)
    }

  }

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

enpls.fit.core = function(plsdf, maxcomp) {

  plsr.cvfit = plsr(y ~ ., data = plsdf,
                    ncomp  = maxcomp,
                    scale  = TRUE,
                    method = 'simpls',
                    validation = 'CV', segments = 5L)

  # select best component number using adjusted CV
  cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])

  plsr.fit = plsr(y ~ ., data = plsdf,
                  ncomp  = cv.bestcomp,
                  scale  = TRUE,
                  method = 'simpls',
                  validation = 'none')

  # save cv.bestcomp for predict.enpls
  enpls.core.fit = list('plsr.fit' = plsr.fit, 'cv.bestcomp' = cv.bestcomp)
  return(enpls.core.fit)

}
