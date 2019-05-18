#' Root Mean Squared Error (RMSE)
#'
#' Compute Root Mean Squared Error (RMSE).
#'
#' @param yreal true response vector
#' @param ypred predicted response vector
#'
#' @return RMSE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export enpls.rmse
enpls.rmse <- function(yreal, ypred) sqrt(mean((yreal - ypred)^2))

#' Mean Absolute Error (MAE)
#'
#' @param yreal true response vector
#' @param ypred predicted response vector
#'
#' @return MAE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export enpls.mae
enpls.mae <- function(yreal, ypred) mean(abs(yreal - ypred))

#' Root Mean Squared Logarithmic Error (RMSLE)
#'
#' @param yreal true response vector
#' @param ypred predicted response vector
#'
#' @return RMSLE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export enpls.rmsle
enpls.rmsle <- function(yreal, ypred)
  sqrt(mean((log(ypred + 1) - log(yreal + 1))^2))
