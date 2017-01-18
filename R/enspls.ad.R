#' Ensemble Sparse Partial Least Squares for
#' Model Applicability Domain Evaluation
#'
#' Model applicability domain evaluation with
#' ensemble sparse partial least squares.
#'
#' @param x Predictor matrix of the training set.
#' @param y Response vector of the training set.
#' @param xtest List, with the i-th component being the i-th test set's
#' predictor matrix (see example code below).
#' @param ytest List, with the i-th component being the i-th test set's
#' response vector (see example code below).
#' @param maxcomp Maximum number of components included within each model.
#' If not specified, will use \code{5} by default.
#' @param cvfolds Number of cross-validation folds used in each model
#' for automatic parameter selection, default is \code{5}.
#' @param alpha Parameter (grid) controlling sparsity of the model.
#' If not specified, default is \code{seq(0.2, 0.8, 0.2)}.
#' @param space Space in which to apply the resampling method.
#' Can be the sample space (\code{"sample"}) or
#' the variable space (\code{"variable"}).
#' @param method Resampling method. \code{"mc"} (Monte-Carlo resampling)
#' or \code{"boot"} (bootstrapping). Default is \code{"mc"}.
#' @param reptimes Number of models to build with Monte-Carlo resampling
#' or bootstrapping.
#' @param ratio Sampling ratio used when \code{method = "mc"}.
#' @param parallel Integer. Number of CPU cores to use.
#' Default is \code{1} (not parallelized).
#'
#' @note Note that for \code{space = "variable"}, \code{method} could
#' only be \code{"mc"}, since bootstrapping in the variable space
#' will create duplicated variables, and that could cause problems.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{tr.error.mean} -
#' absolute mean prediction error for training set
#' \item \code{tr.error.median} -
#' absolute median prediction error for training set
#' \item \code{tr.error.sd} -
#' prediction error sd for training set
#' \item \code{tr.error.matrix} -
#' raw prediction error matrix for training set
#' \item \code{te.error.mean} -
#' list of absolute mean prediction error for test set(s)
#' \item \code{te.error.median} -
#' list of absolute median prediction error for test set(s)
#' \item \code{te.error.sd} -
#' list of prediction error sd for test set(s)
#' \item \code{te.error.matrix} -
#' list of raw prediction error matrix for test set(s)
#' }
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export enspls.ad
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @examples
#' data("logd1k")
#' # remove low variance variables
#' x = logd1k$x[, -c(17, 52, 59)]
#' y = logd1k$y
#'
#' # training set
#' x.tr = x[1:300, ]
#' y.tr = y[1:300]
#'
#' # two test sets
#' x.te = list("test.1" = x[301:400, ],
#'             "test.2" = x[401:500, ])
#' y.te = list("test.1" = y[301:400],
#'             "test.2" = y[401:500])
#'
#' set.seed(42)
#' ad = enspls.ad(x.tr, y.tr, x.te, y.te,
#'                maxcomp = 3, alpha = c(0.3, 0.6, 0.9),
#'                space = "variable", method = "mc",
#'                ratio = 0.8, reptimes = 10)
#' print(ad)
#' plot(ad)
#' # The interactive plot requires a HTML viewer
#' \dontrun{
#' plot(ad, type = "interactive")}

enspls.ad = function(x, y,
                     xtest, ytest,
                     maxcomp  = 5L,
                     cvfolds  = 5L,
                     alpha    = seq(0.2, 0.8, 0.2),
                     space    = c('sample', 'variable'),
                     method   = c('mc', 'boot'),
                     reptimes = 500L,
                     ratio    = 0.8,
                     parallel = 1L) {

  if (missing(x) | missing(y) | missing(xtest) | missing(ytest))
    stop('Please specify x, y, xtest, and ytest')

  if (!all(sapply(xtest, ncol) == ncol(x)))
    stop('All test sets must have identical # of variables as the training set')

  space  = match.arg(space)
  method = match.arg(method)

  n.testset = length(xtest)
  nsamp.tr = nrow(x)
  nsamp.te = sapply(xtest, nrow)

  errorlist.tr = vector('list', reptimes)
  errorlist.te = vector('list', n.testset)
  for (i in 1L:n.testset) errorlist.te[[i]] = vector('list', reptimes)

  if (space == 'sample') {

    idx.row = vector('list', reptimes)

    if (method == 'boot') {

      # 1. space = "sample" & method = "boot"
      for (i in 1L:reptimes) idx.row[[i]] =
          sample(1L:nsamp.tr, nsamp.tr, replace = TRUE)

    } else {

      # 2. space = "sample" & method = "mc"
      n.samp = round(nsamp.tr * ratio)
      for (i in 1L:reptimes) idx.row[[i]] =
          sample(1L:nsamp.tr, n.samp, replace = FALSE)

    }

    if (parallel < 1.5) {

      for (i in 1L:reptimes) {
        fit = enspls.ad.core.fit(x[idx.row[[i]], ], y[idx.row[[i]]],
                                 maxcomp, cvfolds, alpha)
        errorlist.tr[[i]] = enspls.ad.core.pred(fit, x, y)
        for (j in 1L:n.testset) {
          errorlist.te[[j]][[i]] =
            enspls.ad.core.pred(fit, xtest[[j]], ytest[[j]])
        }
      }

    } else {

      registerDoParallel(parallel)
      fit.list = foreach(i = 1L:reptimes) %dopar% {
        enspls.ad.core.fit(x[idx.row[[i]], ], y[idx.row[[i]]],
                           maxcomp, cvfolds, alpha)
      }

      for (i in 1L:reptimes) {
        errorlist.tr[[i]] = enspls.ad.core.pred(fit.list[[i]], x, y)
        for (j in 1L:n.testset) {
          errorlist.te[[j]][[i]] =
            enspls.ad.core.pred(fit.list[[i]], xtest[[j]], ytest[[j]])
        }
      }

    }

  } else {

    if (method == 'boot') {

      # 3. space = "variable" & method = "boot"
      stop('method cannot be "boot" when space = "variable"')

    } else {

      # 4. space = "variable" & method = "mc"
      x.col = ncol(x)
      idx.col = vector('list', reptimes)

      n.var = round(x.col * ratio)
      for (i in 1L:reptimes) idx.col[[i]] =
        sample(1L:x.col, n.var, replace = FALSE)

      if (parallel < 1.5) {

        for (i in 1L:reptimes) {
          fit = enspls.ad.core.fit(x[, idx.col[[i]]], y,
                                   maxcomp, cvfolds, alpha)
          errorlist.tr[[i]] = enspls.ad.core.pred(fit, x[, idx.col[[i]]], y)
          for (j in 1L:n.testset) {
            errorlist.te[[j]][[i]] =
              enspls.ad.core.pred(fit, xtest[[j]][, idx.col[[i]]], ytest[[j]])
          }
        }

      } else {

        registerDoParallel(parallel)
        fit.list = foreach(i = 1L:reptimes) %dopar% {
          enspls.ad.core.fit(x[, idx.col[[i]]], y, maxcomp, cvfolds, alpha)
        }

        for (i in 1L:reptimes) {
          errorlist.tr[[i]] = enspls.ad.core.pred(fit.list[[i]], x[, idx.col[[i]]], y)
          for (j in 1L:n.testset) {
            errorlist.te[[j]][[i]] =
              enspls.ad.core.pred(fit.list[[i]], xtest[[j]][, idx.col[[i]]], ytest[[j]])
          }
        }

      }

    }

  }

  errormat.tr = matrix(NA, ncol = nsamp.tr, nrow = reptimes)
  errormat.te = vector('list', n.testset)
  for (i in 1L:n.testset) errormat.te[[i]] =
    matrix(NA, ncol = nsamp.te[i], nrow = reptimes)

  for (i in 1L:reptimes) errormat.tr[i, ] = as.vector(errorlist.tr[[i]])
  for (i in 1L:reptimes) {
    for (j in 1L:n.testset) {
      errormat.te[[j]][i, ] = as.vector(errorlist.te[[j]][[i]])
    }
  }

  tiny.mean = function (x) abs(colMeans(x, na.rm = TRUE))
  tiny.median = function(x) abs(apply(x, 2L, median, na.rm = TRUE))
  tiny.sd = function(x) apply(x, 2L, sd, na.rm = TRUE)

  tr.error.mean   = tiny.mean(errormat.tr)
  tr.error.median = tiny.mean(errormat.tr)
  tr.error.sd     = tiny.sd(errormat.tr)

  te.error.mean   = lapply(errormat.te, tiny.mean)
  te.error.median = lapply(errormat.te, tiny.median)
  te.error.sd     = lapply(errormat.te, tiny.sd)

  object = list('tr.error.mean'    = tr.error.mean,
                'tr.error.median'  = tr.error.median,
                'tr.error.sd'      = tr.error.sd,
                'tr.error.matrix'  = errormat.tr,
                'te.error.mean'    = te.error.mean,
                'te.error.median'  = te.error.median,
                'te.error.sd'      = te.error.sd,
                'te.error.matrix'  = errormat.te)

  class(object) = 'enspls.ad'
  return(object)

}

#' core fitting function for enspls.ad
#'
#' select the best ncomp and alpha with cross-validation,
#' then use them to fit the complete training set. scale = TRUE
#'
#' @importFrom spls cv.spls spls
#'
#' @return the fitted spls object
#'
#' @keywords internal

enspls.ad.core.fit = function(x.tr, y.tr, maxcomp, cvfolds, alpha) {

  invisible(
    capture.output(
      spls.cvfit <- cv.spls(x.tr,
                            y.tr,
                            fold    = cvfolds,
                            K       = maxcomp,
                            eta     = alpha,
                            scale.x = TRUE,
                            scale.y = FALSE,
                            plot.it = FALSE)))

  # select best component number and alpha using adjusted CV
  cv.bestcomp  = spls.cvfit$'K.opt'
  cv.bestalpha = spls.cvfit$'eta.opt'

  spls.fit = spls(x.tr,
                  y.tr,
                  K       = cv.bestcomp,
                  eta     = cv.bestalpha,
                  scale.x = TRUE,
                  scale.y = FALSE)

  return(spls.fit)

}

#' core prediction function for enspls.ad
#'
#' @return the error vector between predicted and real response
#'
#' @importFrom spls cv.spls spls
#'
#' @keywords internal

enspls.ad.core.pred = function(model, x.te, y.te) {

  pred = predict(model, newx = x.te)

  errorvec = y.te - pred
  names(errorvec) = NULL
  return(errorvec)

}
