#' Make Predictions from a Fitted Ensemble Partial Least Squares Model
#'
#' Make predictions on new data by fitted enpls.fit object.
#'
#' @param object An object of class \code{enpls.fit}.
#' @param newx New data to predict with.
#' @param method Use \code{"mean"} or \code{"median"} to create
#' the final prediction.
#' @param ... Additional parameters for \code{\link{predict}}.
#'
#' @return A numeric vector containing the predicted values.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{enpls.fit}} for fitting ensemble
#' partial least squares regression models.
#'
#' @method predict enpls.fit
#'
#' @export
#'
#' @examples
#' data("alkanes")
#' x <- alkanes$x
#' y <- alkanes$y
#'
#' set.seed(42)
#' fit <- enpls.fit(x, y, reptimes = 50)
#' y.pred <- predict(fit, newx = x)
#' plot(y, y.pred, xlim = range(y), ylim = range(y))
#' abline(a = 0L, b = 1L)
#' y.pred.med <- predict(fit, newx = x, method = "median")
#' plot(y, y.pred.med, xlim = range(y), ylim = range(y))
#' abline(a = 0L, b = 1L)
predict.enpls.fit <- function(object, newx, method = c("mean", "median"), ...) {
  if (missing(newx)) stop("Must provide newx")

  if (!inherits(object, "enpls.fit")) {
    stop('This function only works for objects of class "enpls.fit"')
  }

  method <- match.arg(method)

  nmodel <- length(object)

  predmat <- matrix(NA, ncol = nmodel, nrow = nrow(newx))
  for (i in 1:nmodel) {
    predmat[, i] <- predict(object[[i]][[1]], newx, object[[i]][[2]])
  }

  if (method == "mean") {
    pred <- rowMeans(predmat)
  } else if (method == "median") {
    pred <- apply(predmat, 1L, median)
  }

  pred
}
