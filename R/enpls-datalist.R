#' logD7.4 Data for 1,000 Compounds
#'
#' Distribution coefficients at pH 7.4 (logD7.4) dataset from Wang et, al.
#'
#' This dataset contains distribution coefficients at pH 7.4 (logD7.4)
#' for 1,000 compounds, and 80 molecular descriptors computed with RDKit.
#'
#' @docType data
#' @name logd1k
#' @usage data(logd1k)
#'
#' @format
#' A list with 2 components:
#' \itemize{
#' \item x - data frame with 1,000 rows (samples) and 80 columns (predictors)
#' \item y - numeric vector of length 1,000 (response)
#' }
#' The first 1000 compounds in the original dataset were selected.
#'
#' @references
#' Jian-Bing Wang, Dong-Sheng Cao, Min-Feng Zhu, Yong-Huan Yun, Nan Xiao,
#' and Yi-Zeng Liang. "In silico evaluation of logD7.4 and comparison with
#' other prediction methods."
#' \emph{Journal of Chemometrics} 29, no. 7 (2015): 389--398.
#'
#' @examples
#' data(logd1k)
#' str(logd1k)
NULL

#' Methylalkanes Retention Index Dataset
#'
#' Methylalkanes retention index dataset from Liang et, al.
#'
#' This dataset contains 207 methylalkanes' chromatographic retention index (y)
#' which have been modeled by 21 molecular descriptors (x).
#'
#' Molecular descriptor types:
#' \itemize{
#' \item Chi path, cluster and path/cluster indices
#' \item Kappa shape indices
#' \item E-state indices
#' \item Molecular electricity distance vector index
#' }
#'
#' @docType data
#' @name alkanes
#' @usage data("alkanes")
#'
#' @format
#' A list with 2 components:
#' \itemize{
#' \item x - data frame with 207 rows (samples) and 21 columns (predictors)
#' \item y - numeric vector of length 207 (response)
#' }
#'
#' @references
#' Yi-Zeng Liang, Da-Lin Yuan, Qing-Song Xu, and Olav Martin Kvalheim.
#' "Modeling based on subspace orthogonal projections for QSAR and
#' QSPR research." \emph{Journal of Chemometrics} 22, no. 1 (2008): 23--35.
#'
#' @examples
#' data("alkanes")
#' str(alkanes)
NULL
