#' Make Predictions from a Fitted Sparse Ensemble Partial Least Squares Model
#'
#' Make predictions on new data by fitted enspls.fit object.
#'
#' @param object An object of class \code{enspls.fit}.
#' @param newx New data to predict with.
#' @param method Use \code{mean} or \code{median} as the final prediction.
#' @param ... Additional parameters for \code{\link{predict}}.
#'
#' @return A numeric vector containing the predicted values.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{enspls.fit}} for ensemble sparse
#' partial least squares regression.
#'
#' @method predict enspls.fit
#'
#' @importFrom spls predict.spls
#'
#' @export
#'
#' @examples
#' data("logd1k")
#' x = logd1k$x
#' y = logd1k$y
#'
#' set.seed(42)
#' fit = enspls.fit(x, y, MCtimes = 5, maxcomp = 2)
#' y.pred = predict(fit, newx = x)
#' plot(y, y.pred, xlim = range(y), ylim = range(y))
#' abline(a = 0L, b = 1L)
#' y.pred.med = predict(fit, newx = x, method = 'median')
#' plot(y, y.pred.med, xlim = range(y), ylim = range(y))
#' abline(a = 0L, b = 1L)

predict.enspls.fit = function(object, newx, method = c('mean', 'median'), ...) {

  if (missing(newx)) stop('Must provide newx')

  if (!inherits(object, 'enspls.fit'))
    stop('This function only works for objects of class "enspls.fit"')

  method = match.arg(method)

  nmodel = length(object)

  predmat = matrix(NA, ncol = nmodel, nrow = nrow(newx))
  for (i in 1:nmodel) {
    predmat[, i] = predict(object[[i]][[1]], newx)
  }

  if (method == 'mean') {
    pred = rowMeans(predmat)
  } else if (method == 'median') {
    pred = apply(predmat, 1L, median)
  }

  return(pred)

}
