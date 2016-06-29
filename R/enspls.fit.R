#' Ensemble Sparse Partial Least Squares Regression
#'
#' Ensemble sparse partial least squares regression.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param maxcomp Maximum number of components included within the models,
#' if not specified, default is 5.
#' @param alpha Parameter (grid) controlling sparsity of the model.
#' If not specified, default is \code{seq(0.2, 0.8, 0.2)}.
#' @param MCtimes times of Monte-Carlo
#' @param method \code{"mc"} or \code{"bootstrap"}. Default is \code{"mc"}.
#' @param ratio sample ratio used when \code{method = "mc"}
#' @param parallel Integer. Number of CPU cores to use.
#' Default is \code{1} (not parallelized).
#'
#' @return A list containing all sparse partial least squares model objects.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enspls.fs}} for feature selection with ensemble
#' sparse partial least squares regression.
#' See \code{\link{enspls.od}} for outlier detection with ensemble
#' sparse partial least squares regression.
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
#' fit = enspls.fit(x, y, MCtimes = 4, maxcomp = 3)
#' print(fit)
#' predict(fit, newx = x)

enspls.fit = function(x, y,
                      maxcomp = 5L,
                      alpha = seq(0.2, 0.8, 0.2),
                      MCtimes = 500L,
                      method = c('mc', 'bootstrap'), ratio = 0.8,
                      parallel = 1L) {

  if (missing(x) | missing(y)) stop('Please specify both x and y')

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
      modellist[[i]] = enspls.fit.core(xtmp, ytmp, maxcomp, alpha)
    }

  } else {

    registerDoParallel(parallel)
    modellist = foreach(i = 1L:MCtimes) %dopar% {
      xtmp = x[samp.idx[[i]], ]
      ytmp = y[samp.idx[[i]]]
      enspls.fit.core(xtmp, ytmp, maxcomp, alpha)
    }

  }

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

enspls.fit.core = function(xtmp, ytmp, maxcomp, alpha) {

  invisible(
    capture.output(
      spls.cvfit <- cv.spls(xtmp, ytmp,
                            fold = 5,
                            K = maxcomp, eta = alpha,
                            scale.x = TRUE, scale.y = FALSE,
                            plot.it = FALSE)))

  # select best component number and alpha using adjusted CV
  cv.bestcomp  = spls.cvfit$'K.opt'
  cv.bestalpha = spls.cvfit$'eta.opt'

  spls.fit = spls(xtmp, ytmp,
                  K = cv.bestcomp, eta = cv.bestalpha,
                  scale.x = TRUE, scale.y = FALSE)

  # save cv.bestcomp and cv.bestalpha for predict.enspls
  enspls.core.fit = list('spls.fit' = spls.fit,
                         'cv.bestcomp' = cv.bestcomp,
                         'cv.bestalpha' = cv.bestalpha)
  return(enspls.core.fit)

}
