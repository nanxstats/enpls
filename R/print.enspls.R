#' Print cv.enspls Object
#'
#' Print cv.enspls object.
#'
#' @param x An object of class \code{cv.enspls}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{cv.enspls}} for cross-validation of
#' ensemble sparse partial least squares regression models.
#'
#' @method print cv.enspls
#'
#' @export
#'
#' @examples
#' # This example takes one minute to run
#' \dontrun{
#' data("logd1k")
#' x = logd1k$x
#' y = logd1k$y
#'
#' set.seed(42)
#' cvfit = cv.enspls(x, y, reptimes = 10)
#' print(cvfit)}

print.cv.enspls = function(x, ...) {

  if (!inherits(x, 'cv.enspls'))
    stop('This function only works for objects of class "cv.enspls"')

  cat('Cross Validation Result for Ensemble Sparse Partial Least Squares\n')
  cat('---\n')
  cat(paste(
    'RMSE = ', sprintf("%.4f", x$'RMSE'), '\n',
    'MAE = ', sprintf("%.6f", x$'MAE'), '\n',
    'Rsquare = ', sprintf("%.6f", x$'Rsquare'), '\n',
    sep = ''))

}

#' Print Fitted Ensemble Sparse Partial Least Squares Object
#'
#' Print coefficients of each model in the enspls.fit object.
#'
#' @param x An object of class \code{enspls.fit}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{enspls.fit}} for fitting ensemble sparse
#' partial least squares regression models.
#'
#' @method print enspls.fit
#'
#' @importFrom spls coef.spls
#'
#' @export
#'
#' @examples
#' data("logd1k")
#' x = logd1k$x
#' y = logd1k$y
#'
#' set.seed(42)
#' fit = enspls.fit(
#'   x, y, reptimes = 5, maxcomp = 3,
#'   alpha = c(0.3, 0.6, 0.9))
#' print(fit)

print.enspls.fit = function(x, ...) {

  if (!inherits(x, 'enspls.fit'))
    stop('This function only works for objects of class "enspls.fit"')

  coefmeta = coef(x[[1]][[1]])
  varcount = nrow(coefmeta)
  reptimes  = length(x)
  coefdf   = matrix(NA, ncol = reptimes, nrow = varcount)
  for (i in 1:reptimes) coefdf[, i] = coef(x[[i]][[1]])
  rownames(coefdf) = rownames(coefmeta)

  cat('Coefficients of the Models by Ensemble Sparse Partial Least Squares\n')
  cat('---\n')
  print(coefdf)

}

#' Print enspls.fs Object
#'
#' Print enspls.fs object.
#'
#' @param x An object of class \code{enspls.fs}.
#' @param sort Should the variables be sorted in decreasing order of importance?
#' @param nvar Number of top variables to show. Ignored if \code{sort = FALSE}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{enspls.fs}} for measuring feature importance with
#' ensemble sparse partial least squares regressions.
#'
#' @method print enspls.fs
#'
#' @export
#'
#' @examples
#' data("logd1k")
#' x = logd1k$x
#' y = logd1k$y
#'
#' set.seed(42)
#' fs = enspls.fs(
#'   x, y, reptimes = 5, maxcomp = 3,
#'   alpha = c(0.3, 0.6, 0.9))
#' print(fs, nvar = 10L)

print.enspls.fs = function(x, sort = TRUE, nvar = NULL, ...) {

  if (!inherits(x, 'enspls.fs'))
    stop('This function only works for objects of class "enspls.fs"')

  varimp = x$'variable.importance'
  if (is.null(nvar)) nvar = length(varimp)

  cat('Variable Importance by Ensemble Sparse Partial Least Squares\n')
  cat('---\n')
  if (sort == TRUE) {
    print(data.frame('Importance' = sort(varimp, TRUE)[1:nvar]))
  } else {
    print(data.frame('Importance' = varimp))
  }

}

#' Print enspls.od Object
#'
#' Print enspls.od object.
#'
#' @param x An object of class \code{enspls.od}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{enspls.od}} for outlier detection with
#' ensemble sparse partial least squares regressions.
#'
#' @method print enspls.od
#'
#' @export
#'
#' @examples
#' data("logd1k")
#' x = logd1k$x
#' y = logd1k$y
#'
#' set.seed(42)
#' od = enspls.od(
#'   x, y, reptimes = 5, maxcomp = 3,
#'   alpha = c(0.3, 0.6, 0.9))
#' print(od)

print.enspls.od = function(x, ...) {

  if (!inherits(x, 'enspls.od'))
    stop('This function only works for objects of class "enspls.od"')

  cat('Outlier Detection by Ensemble Sparse Partial Least Squares\n')
  cat('---\n')
  cat('Mean residual for each sample:\n')
  print(x$'error.mean')
  cat('---\n')
  cat('Residual SD for each sample:\n')
  print(x$'error.sd')

}

#' Print enspls.ad Object
#'
#' Print enspls.ad object.
#'
#' @param x An object of class \code{enspls.ad}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{enspls.ad}} for model applicability domain
#' evaluation with ensemble sparse partial least squares regressions.
#'
#' @method print enspls.ad
#'
#' @export
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
#' x.te = list(
#'   "test.1" = x[301:400, ],
#'   "test.2" = x[401:500, ])
#' y.te = list(
#'   "test.1" = y[301:400],
#'   "test.2" = y[401:500])
#'
#' set.seed(42)
#' ad = enspls.ad(
#'   x.tr, y.tr, x.te, y.te,
#'   maxcomp = 3, alpha = c(0.3, 0.6, 0.9),
#'   space = "variable", method = "mc",
#'   ratio = 0.8, reptimes = 10)
#' print(ad)

print.enspls.ad = function(x, ...) {

  if (!inherits(x, 'enspls.ad'))
    stop('This function only works for objects of class "enspls.ad"')

  cat('Model Applicability Domain Evaluation by ENSPLS\n')
  cat('---\n')
  cat('Absolute mean prediction error for each training set sample:\n')
  print(x$'tr.error.mean')
  cat('---\n')
  cat('Prediction error SD for each training set sample:\n')
  print(x$'tr.error.sd')
  cat('---\n')
  cat('Absolute mean prediction error for each test set sample:\n')
  print(x$'te.error.mean')
  cat('---\n')
  cat('Prediction error SD for each test set sample:\n')
  print(x$'te.error.sd')

}
