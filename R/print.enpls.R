#' Print cv.enpls Object
#'
#' Print cv.enpls object.
#'
#' @param x An object of class \code{cv.enpls}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{cv.enpls}} for ensemble
#' partial least squares regression.
#'
#' @method print cv.enpls
#'
#' @export
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' cvfit = cv.enpls(x, y, MCtimes = 10)
#' print(cvfit)

print.cv.enpls = function(x, ...) {

  if (!inherits(x, 'cv.enpls'))
    stop('This function only works for objects of class "cv.enpls"')

  cat('Cross Validation Result for Ensemble Partial Least Squares\n')
  cat('---\n')
  cat(paste('RMSE = ', sprintf("%.4f", x$'RMSE'), '\n',
            'MAE = ', sprintf("%.6f", x$'MAE'), '\n',
            'Rsquare = ', sprintf("%.6f", x$'Rsquare'), '\n',
            sep = ''))

}

#' Print Fitted Ensemble Partial Least Squares Object
#'
#' Print coefficients of each model in the enpls.fit object.
#'
#' @param x An object of class \code{enpls.fit}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enpls.fit}} for ensemble
#' partial least squares regression.
#'
#' @method print enpls.fit
#'
#' @export
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' fit = enpls.fit(x, y, MCtimes = 50)
#' print(fit)

print.enpls.fit = function(x, ...) {

  if (!inherits(x, 'enpls.fit'))
    stop('This function only works for objects of class "enpls.fit"')

  coefmeta = coef(x[[1]][[1]], intercept = TRUE)[, 1, 1]
  varcount = length(coefmeta)
  mctimes  = length(x)
  coefdf   = matrix(NA, ncol = mctimes, nrow = varcount)
  for (i in 1:mctimes) coefdf[, i] = coef(x[[i]][[1]], intercept = TRUE)[, 1, 1]
  rownames(coefdf) = names(coefmeta)

  cat('Coefficients of the Models by Ensemble Partial Least Squares\n')
  cat('---\n')
  print(coefdf)

}

#' Print enpls.fs Object
#'
#' Print enpls.fs object.
#'
#' @param x An object of class \code{enpls.fs}.
#' @param sort Should the variables be sorted in decreasing order of importance?
#' @param nvar How many variables to show? Ignored if \code{sort = FALSE}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enpls.fs}} for feature selection with
#' ensemble partial least squares regression.
#'
#' @method print enpls.fs
#'
#' @export
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' fs = enpls.fs(x, y, MCtimes = 100)
#' print(fs)
#' print(fs, nvar = 10L)

print.enpls.fs = function(x, sort = TRUE, nvar = NULL, ...) {

  if (!inherits(x, 'enpls.fs'))
    stop('This function only works for objects of class "enpls.fs"')

  varimp = x$'variable.importance'
  if (is.null(nvar)) nvar = length(varimp)

  cat('Variable Importance by Ensemble Partial Least Squares\n')
  cat('---\n')
  if (sort == TRUE) {
    print(data.frame('Importance' = sort(varimp, TRUE)[1:nvar]))
  } else {
    print(data.frame('Importance' = varimp))
  }

}

#' Print enpls.od Object
#'
#' Print enpls.od object.
#'
#' @param x An object of class \code{enpls.od}.
#' @param ... Additional parameters for \code{\link{print}}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enpls.od}} for outlier detection with
#' ensemble partial least squares regression.
#'
#' @method print enpls.od
#'
#' @export
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' od = enpls.od(x, y, MCtimes = 50)
#' print(od)

print.enpls.od = function(x, ...) {

  if (!inherits(x, 'enpls.od'))
    stop('This function only works for objects of class "enpls.od"')

  cat('Outlier Detection by Ensemble Partial Least Squares\n')
  cat('---\n')
  cat('Mean residual for each sample:\n')
  print(x$'error.mean')
  cat('---\n')
  cat('Residual SD for each sample:\n')
  print(x$'error.sd')

}
