#' Cross Validation for Ensemble Partial Least Squares Regression
#'
#' K-fold cross validation for ensemble partial least squares regression.
#'
#' @param x Predictor matrix.
#' @param y Response vector.
#' @param nfolds Number of cross-validation folds, default is \code{5}.
#' Note that this is the CV folds for the ensemble PLS model,
#' not the individual PLS models. To control the CV folds for
#' single PLS models, please use the argument \code{cvfolds}.
#' @param verbose Shall we print out the progress of cross-validation?
#' @param ... Arguments to be passed to \code{\link{enpls.fit}}.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{ypred} - a matrix containing two columns: real y and predicted y
#' \item \code{residual} - cross validation result (y.pred - y.real)
#' \item \code{RMSE} - RMSE
#' \item \code{MAE} - MAE
#' \item \code{Rsquare} - Rsquare
#' }
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @note To maximize the probablity that each observation can
#' be selected in the test set (thus the prediction uncertainty
#' can be measured), please try setting a large \code{reptimes}.
#'
#' @seealso See \code{\link{enpls.fit}} for ensemble
#' partial least squares regressions.
#'
#' @export cv.enpls
#'
#' @examples
#' data("alkanes")
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' cvfit = cv.enpls(x, y, reptimes = 10)
#' print(cvfit)
#' plot(cvfit)

cv.enpls = function(x, y, nfolds = 5L, verbose = TRUE, ...) {

  if (missing(x) | missing(y)) stop('Please specify both x and y')

  x.row = nrow(x)
  index = rep_len(1L:nfolds, x.row)

  ypred = matrix(NA, ncol = 2L, nrow = x.row)

  for (i in 1L:nfolds) {
    if (verbose) cat('Beginning fold', i, '\n')
    xtrain = x[index != i, ]
    ytrain = y[index != i]
    xtest  = x[index == i, ]
    ytest  = y[index == i]
    fit = enpls.fit(xtrain, ytrain, ...)
    ypredvec = predict(fit, newx = xtest)
    ypred[index == i, 1L] = ytest
    ypred[index == i, 2L] = ypredvec
  }

  colnames(ypred) = c('y.real', 'y.pred')

  residual = ypred[, 1L] - ypred[, 2L]
  RMSE     = sqrt(mean((residual)^2, na.rm = TRUE))
  MAE      = mean(abs(residual), na.rm = TRUE)
  Rsquare  = 1L - (sum((residual)^2, na.rm = TRUE)/sum((y - mean(y))^2))

  object = list('ypred'    = ypred,
                'residual' = residual,
                'RMSE'     = RMSE,
                'MAE'      = MAE,
                'Rsquare'  = Rsquare)
  class(object) = 'cv.enpls'
  return(object)

}
