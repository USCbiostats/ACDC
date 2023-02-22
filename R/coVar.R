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
#' @export
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