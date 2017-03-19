#' Ensemble Partial Least Squares for Model Applicability Domain Evaluation
#'
#' Model applicability domain evaluation with ensemble partial least squares.
#'
#' @param x Predictor matrix of the training set.
#' @param y Response vector of the training set.
#' @param xtest List, with the i-th component being the i-th test set's
#' predictor matrix (see example code below).
#' @param ytest List, with the i-th component being the i-th test set's
#' response vector (see example code below).
#' @param maxcomp Maximum number of components included within each model.
#' If not specified, will use the maximum number possible (considering
#' cross-validation and special cases where n is smaller than p).
#' @param cvfolds Number of cross-validation folds used in each model
#' for automatic parameter selection, default is \code{5}.
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
#' @export enpls.ad
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' # training set
#' x.tr = x[1:100, ]
#' y.tr = y[1:100]
#'
#' # two test sets
#' x.te = list(
#'   "test.1" = x[101:150, ],
#'   "test.2" = x[151:207, ])
#' y.te = list(
#'   "test.1" = y[101:150],
#'   "test.2" = y[151:207])
#'
#' set.seed(42)
#' ad = enpls.ad(
#'   x.tr, y.tr, x.te, y.te,
#'   space = "variable", method = "mc",
#'   ratio = 0.9, reptimes = 50)
#' print(ad)
#' plot(ad)
#' # the interactive plot requires a HTML viewer
#' \dontrun{
#' plot(ad, type = "interactive")}

enpls.ad = function(
  x, y,
  xtest, ytest,
  maxcomp  = NULL,
  cvfolds  = 5L,
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
        plsdf.tr = as.data.frame(cbind(x[idx.row[[i]], ],
                                       'y' = y[idx.row[[i]]]))
        fit = suppressWarnings(enpls.ad.core.fit(plsdf.tr, maxcomp, cvfolds))
        errorlist.tr[[i]] = suppressWarnings(enpls.ad.core.pred(
          fit, as.data.frame(cbind(x, 'y' = y))))
        for (j in 1L:n.testset) {
          errorlist.te[[j]][[i]] =
            suppressWarnings(enpls.ad.core.pred(
              fit, as.data.frame(cbind(xtest[[j]], 'y' = ytest[[j]]))))
        }
      }

    } else {

      registerDoParallel(parallel)
      fit.list = foreach(i = 1L:reptimes) %dopar% {
        plsdf.tr = as.data.frame(cbind(x[idx.row[[i]], ], 'y' = y[idx.row[[i]]]))
        enpls.ad.core.fit(plsdf.tr, maxcomp, cvfolds)
      }

      for (i in 1L:reptimes) {
        errorlist.tr[[i]] = suppressWarnings(enpls.ad.core.pred(
          fit.list[[i]], as.data.frame(cbind(x, 'y' = y))))
        for (j in 1L:n.testset) {
          errorlist.te[[j]][[i]] =
            suppressWarnings(enpls.ad.core.pred(
              fit.list[[i]], as.data.frame(cbind(xtest[[j]], 'y' = ytest[[j]]))))
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
          plsdf.tr = as.data.frame(cbind(x[, idx.col[[i]]], 'y' = y))
          fit = suppressWarnings(enpls.ad.core.fit(plsdf.tr, maxcomp, cvfolds))
          errorlist.tr[[i]] = suppressWarnings(enpls.ad.core.pred(fit, plsdf.tr))
          for (j in 1L:n.testset) {
            errorlist.te[[j]][[i]] =
              suppressWarnings(enpls.ad.core.pred(
                fit, as.data.frame(cbind(xtest[[j]][, idx.col[[i]]], 'y' = ytest[[j]]))))
          }
        }

      } else {

        registerDoParallel(parallel)
        fit.list = foreach(i = 1L:reptimes) %dopar% {
          plsdf.tr = as.data.frame(cbind(x[, idx.col[[i]]], 'y' = y))
          enpls.ad.core.fit(plsdf.tr, maxcomp, cvfolds)
        }

        for (i in 1L:reptimes) {
          errorlist.tr[[i]] = suppressWarnings(enpls.ad.core.pred(
            fit.list[[i]], as.data.frame(cbind(x[, idx.col[[i]]], 'y' = y))))
          for (j in 1L:n.testset) {
            errorlist.te[[j]][[i]] =
              suppressWarnings(enpls.ad.core.pred(
                fit.list[[i]], as.data.frame(cbind(xtest[[j]][, idx.col[[i]]], 'y' = ytest[[j]]))))
          }
        }

      }

    }

  }

  errormat.tr = matrix(NA, ncol = nsamp.tr, nrow = reptimes)
  errormat.te = vector('list', n.testset)
  for (i in 1L:n.testset) errormat.te[[i]] =
    matrix(NA, ncol = nsamp.te[i], nrow = reptimes)

  for (i in 1L:reptimes) errormat.tr[i, ] = errorlist.tr[[i]]
  for (i in 1L:reptimes) {
    for (j in 1L:n.testset) {
      errormat.te[[j]][i, ] = errorlist.te[[j]][[i]]
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

  res = list(
    'tr.error.mean'    = tr.error.mean,
    'tr.error.median'  = tr.error.median,
    'tr.error.sd'      = tr.error.sd,
    'tr.error.matrix'  = errormat.tr,
    'te.error.mean'    = te.error.mean,
    'te.error.median'  = te.error.median,
    'te.error.sd'      = te.error.sd,
    'te.error.matrix'  = errormat.te)
  class(res) = 'enpls.ad'

  res

}

#' core fitting function for enpls.ad
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set. scale = TRUE
#'
#' @return the fitted plsr object
#'
#' @keywords internal

enpls.ad.core.fit = function(trainingset, maxcomp, cvfolds) {

  if (is.null(maxcomp)) {

    plsr.cvfit = plsr(
      y ~ .,
      data       = trainingset,
      scale      = TRUE,
      method     = 'simpls',
      validation = 'CV',
      segments   = cvfolds)

  } else {

    plsr.cvfit = plsr(
      y ~ .,
      data       = trainingset,
      ncomp      = maxcomp,
      scale      = TRUE,
      method     = 'simpls',
      validation = 'CV',
      segments   = cvfolds)

  }

  # select best component number using adjusted CV
  cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])

  plsr.fit = plsr(
    y ~ .,
    data       = trainingset,
    ncomp      = cv.bestcomp,
    scale      = TRUE,
    method     = 'simpls',
    validation = 'none')

  list('plsr.fit' = plsr.fit, 'cv.bestcomp' = cv.bestcomp)

}

#' core prediction function for enpls.ad
#'
#' @return the error vector between predicted and real response
#'
#' @keywords internal

enpls.ad.core.pred = function(model, testset) {

  pred = predict(
    model$'plsr.fit', ncomp = model$'cv.bestcomp',
    newdata = testset[, !(colnames(testset) %in% c('y'))])[, 1L, 1L]

  errorvec = testset[, 'y'] - pred
  names(errorvec) = NULL

  errorvec

}
