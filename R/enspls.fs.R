#' Ensemble Sparse Partial Least Squares for Feature Selection
#'
#' Feature selection with ensemble sparse partial least squares.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param maxcomp Maximum number of components included within the models,
#' if not specified, default is 5.
#' @param alpha Parameter (grid) controlling sparsity of the model.
#' If not specified, default is \code{seq(0.2, 0.8, 0.2)}.
#' @param reptimes Number of models to build with Monte-Carlo resampling
#' or bootstrapping.
#' @param method \code{"mc"} or \code{"bootstrap"}. Default is \code{"mc"}.
#' @param ratio sample ratio used when \code{method = "mc"}
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
#' @seealso See \code{\link{enspls.od}} for outlier detection with
#' ensemble sparse partial least squares regression.
#' See \code{\link{enspls.fit}} for ensemble sparse
#' partial least squares regression.
#'
#' @export enspls.fs
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
#' fs = enspls.fs(x, y, reptimes = 5, maxcomp = 2)
#' print(fs, nvar = 10)
#' plot(fs, nvar = 10)
#' plot(fs, type = 'boxplot', limits = c(0.05, 0.95), nvar = 10)

enspls.fs = function(x, y,
                     maxcomp = 5L,
                     alpha = seq(0.2, 0.8, 0.2),
                     reptimes = 500L,
                     method = c('mc', 'bootstrap'), ratio = 0.8,
                     parallel = 1L) {

  if (missing(x) | missing(y)) stop('Please specify both x and y')

  method = match.arg(method)

  x.row = nrow(x)
  samp.idx = vector('list', reptimes)

  if (method == 'mc') {
    for (i in 1L:reptimes) samp.idx[[i]] = sample(1L:x.row, round(x.row * ratio))
  }

  if (method == 'bootstrap') {
    for (i in 1L:reptimes) samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
  }

  if (parallel < 1.5) {

    coeflist = vector('list', reptimes)
    for (i in 1L:reptimes) {
      xtmp = x[samp.idx[[i]], ]
      xtmp = scale(xtmp, center = TRUE, scale = TRUE)
      ytmp = y[samp.idx[[i]]]
      coeflist[[i]] = enspls.fs.core(xtmp, ytmp, maxcomp, alpha)
    }

  } else {

    registerDoParallel(parallel)
    coeflist = foreach(i = 1L:reptimes) %dopar% {
      xtmp = x[samp.idx[[i]], ]
      xtmp = scale(xtmp, center = TRUE, scale = TRUE)
      ytmp = y[samp.idx[[i]]]
      enspls.fs.core(xtmp, ytmp, maxcomp, alpha)
    }

  }

  coefmat = do.call(rbind, coeflist)

  varimp = abs(colMeans(coefmat))/apply(coefmat, 2L, sd)

  # let variables with all zero coefficients have 0 importance
  varimp[which(!is.finite(varimp))] = 0

  object = list('variable.importance' = varimp,
                'coefficient.matrix'  = coefmat)
  class(object) = 'enspls.fs'
  return(object)

}

#' core function for enspls.fs
#'
#' select the best ncomp and alpha with cross-validation, then
#' use it to fit the complete training set.
#' scale = FALSE
#'
#' @importFrom spls cv.spls spls
#'
#' @return fitted coefficients
#'
#' @keywords internal

enspls.fs.core = function(xtmp, ytmp, maxcomp, alpha) {

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

  spls.coef = coef(spls.fit)
  spls.coef.vec = as.vector(spls.coef)
  names(spls.coef.vec) = rownames(spls.coef)

  return(spls.coef.vec)

}
