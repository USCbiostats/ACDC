#' coVar
#'
#' @description Function to calculate ACDC covariances within a data pair for all samples
#' 
#' @param dataPair column indices of two genes to calculate covariance between
#' @param fullData dataframe or matrix with samples as rows, all probes as columns
#' @return Vector of pairwise covariances for all samples
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @examples
#' #load Biobase package
#' library(Biobase)
#' 
#' # load the dataset
#' data(sample.ExpressionSet)
#' 
#' # pick 100 random samples
#' samp <- sample(1:500, size=100)
#' 
#' # first 100 probes from expression data, transposed so samples are rows
#' data <- t(exprs(sample.ExpressionSet)[samp,])
#' 
#' # run function with first two samples
#' coVar(dataPair = c(1, 2), 
#'       fullData = data)
#' 
#' @export
#' @import Biobase
coVar <- function(dataPair, fullData) {
  cv             <- numeric(nrow(fullData))
  meanExpression <- apply(fullData, MARGIN = 2, FUN = mean)
  
  for (l in 1:nrow(fullData)) {
    g_j    <- fullData[l, dataPair[1]]
    g_k    <- fullData[l, dataPair[2]]
    bar_gj <- meanExpression[dataPair[1]]
    bar_gk <- meanExpression[dataPair[2]]
    cv[l]  <- ((g_j-bar_gj)*(g_k-bar_gk))
  }
  
  return(cv)
}