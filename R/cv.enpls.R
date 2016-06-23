#' Cross Validation for Ensemble Partial Least Squares Regression
#'
#' K-fold cross validation for ensemble partial least squares regression.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param nfolds number of folds - default is \code{5}.
#' @param verbose shall we print the cross validation process
#' @param ... other arguments that can be passed to \code{\link{enpls.fit}}
#'
#' @return A list containing four components:
#' \itemize{
#' \item \code{ypred} - a matrix containing two columns: real y and predicted y
#' \item \code{residual} - cross validation result (y.pred - y.real)
#' \item \code{RMSE} - RMSE
#' \item \code{R2} - R2
#' }
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @note To maximize the probablity that each observation can
#' be selected in the test set (thus the prediction uncertainty
#' can be measured), please try setting a large \code{MCtimes}.
#'
#' @seealso See \code{\link{enpls.fit}} for ensemble
#' partial least squares regression.
#'
#' @export cv.enpls
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
#' cvfit = cv.enpls(x, y, MCtimes = 20)
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
    ypredvec  = predict(fit, newx = xtest)
    ypred[index == i, 1L] = ytest
    ypred[index == i, 2L] = ypredvec
  }

  colnames(ypred) = c('y.real', 'y.pred')

  residual = ypred[, 1L] - ypred[, 2L]
  RMSE = sqrt(mean((residual)^2, na.rm = TRUE))
  R2 = 1L - (sum((residual)^2, na.rm = TRUE)/sum((y - mean(y))^2))

  object = list('ypred' = ypred,
                'residual' = residual,
                'RMSE' = RMSE,
                'R2' = R2)
  class(object) = 'cv.enpls'
  return(object)

}
