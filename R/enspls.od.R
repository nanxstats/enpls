#' Ensemble Sparse Partial Least Squares for Outlier Detection
#'
#' Outlier detection with ensemble sparse partial least squares.
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
#' @return A list containing four components:
#' \itemize{
#' \item \code{error.mean} - error mean for all samples (absolute value)
#' \item \code{error.median} - error median for all samples
#' \item \code{error.sd} - error sd for all samples
#' \item \code{predict.error.matrix} - the original prediction error matrix
#' }
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @note To maximize the probablity that each observation can
#' be selected in the test set (thus the prediction uncertainty
#' can be measured), please try setting a large \code{reptimes}.
#'
#' @seealso See \code{\link{enspls.fs}} for measuring feature importance
#' with ensemble sparse partial least squares regressions.
#' See \code{\link{enspls.fit}} for fitting ensemble sparse
#' partial least squares regression models.
#'
#' @export enspls.od
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
#' od = enspls.od(x, y, reptimes = 5, maxcomp = 3,
#'                alpha = c(0.3, 0.6, 0.9))
#' plot(od, prob = 0.1)
#' plot(od, criterion = "sd", sdtimes = 1)

enspls.od = function(x, y,
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
  samp.idx.remain = vector('list', reptimes)

  if (method == 'mc') {
    for (i in 1L:reptimes) {
      samp.idx[[i]] = sample(1L:x.row, round(x.row * ratio))
      samp.idx.remain[[i]] = setdiff(1L:x.row, samp.idx[[i]])
    }
  }

  if (method == 'boot') {
    for (i in 1L:reptimes) {
      samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
      samp.idx.remain[[i]] = setdiff(1L:x.row, unique(samp.idx[[i]]))
    }
  }

  if (parallel < 1.5) {

    errorlist = vector('list', reptimes)
    for (i in 1L:reptimes) {
      x.sample = x[samp.idx[[i]], ]
      x.remain = x[samp.idx.remain[[i]], ]
      y.sample = y[samp.idx[[i]]]
      y.remain = y[samp.idx.remain[[i]]]
      errorlist[[i]] =
        enspls.od.core(x.sample, y.sample, x.remain, y.remain,
                       maxcomp, cvfolds, alpha)
    }

  } else {

    registerDoParallel(parallel)
    errorlist = foreach(i = 1L:reptimes) %dopar% {
      x.sample = x[samp.idx[[i]], ]
      x.remain = x[samp.idx.remain[[i]], ]
      y.sample = y[samp.idx[[i]]]
      y.remain = y[samp.idx.remain[[i]]]
      enspls.od.core(x.sample, y.sample, x.remain, y.remain,
                     maxcomp, cvfolds, alpha)
    }

  }

  prederrmat = matrix(NA, ncol = x.row, nrow = reptimes)
  for (i in 1L:reptimes) {
    for (j in 1L:length(samp.idx.remain[[i]])) {
      prederrmat[i, samp.idx.remain[[i]][j]] = errorlist[[i]][j]
    }
  }

  errmean   = abs(colMeans(prederrmat, na.rm = TRUE))
  errmedian = apply(prederrmat, 2L, median, na.rm = TRUE)
  errsd     = apply(prederrmat, 2L, sd, na.rm = TRUE)

  object = list('error.mean'    = errmean,
                'error.median'  = errmedian,
                'error.sd'      = errsd,
                'predict.error.matrix' = prederrmat)
  class(object) = 'enspls.od'
  return(object)

}

#' core function for enspls.od
#'
#' select the best ncomp and alpha with cross-validation,
#' then use them to fit the complete training set,
#' and predict on the test set. scale = TRUE
#'
#' @importFrom spls cv.spls spls
#'
#' @return the error vector between predicted y and real y
#'
#' @keywords internal

enspls.od.core = function(x.sample, y.sample, x.remain, y.remain,
                          maxcomp, cvfolds, alpha) {

  invisible(
    capture.output(
      spls.cvfit <- cv.spls(x.sample,
                            y.sample,
                            fold    = cvfolds,
                            K       = maxcomp,
                            eta     = alpha,
                            scale.x = TRUE,
                            scale.y = FALSE,
                            plot.it = FALSE)))

  # select best component number and alpha using adjusted CV
  cv.bestcomp  = spls.cvfit$'K.opt'
  cv.bestalpha = spls.cvfit$'eta.opt'

  spls.fit = spls(x.sample,
                  y.sample,
                  K       = cv.bestcomp,
                  eta     = cv.bestalpha,
                  scale.x = TRUE,
                  scale.y = FALSE)

  pred = predict(spls.fit, newx = x.remain)

  error = y.remain - pred
  names(error) = NULL

  return(error)

}
