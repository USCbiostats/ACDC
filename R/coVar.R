#' coVar
#'
#' @description Function to calculate ACDC covariances within a data pair for all samples
#' 
#' @param dataPair column indices of two genes to calculate covariance between
#' @param fullData dataframe or matrix with samples as rows, all probes as columns
#' @return Co-expression profile, or pairwise covariances for all samples, vector for given features
#' 
#' @details Co-expression for a single sample, s, is defined as
#' \deqn{c_{s,j,k} \equiv \left(g_{s,j}-\bar{g_j}\right)\left(g_{s,k}-\bar{g_k}\right)}
#' where \eqn{g_{s,j}} denotes the expression of gene j in sample s and \eqn{\bar{g_j}} denotes the mean expression of gene j in all samples.
#' 
#' Denoting the sample size as N, coVar returns the co-expression profile across all samples:
#' \deqn{c_{j,k} = (c_{1,j,k}, c_{2,j,k}, ... , c_{N,j,k})}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @references 
#' Martin P, et al. Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed through a nutrigenomic study. Hepatology, in press, 2007.
#' 
#' @examples
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # run function with first two samples
#' coVar(dataPair = c(1, 2), 
#'       fullData = nutrimouse$lipid)
#' 
#' @export
coVar <- function(dataPair, fullData) {
  
  meanExpression <- apply(fullData, MARGIN = 2, FUN = mean)
  
  # helper function to calculate covariance
  coVar_calc <- function(x) {
    g_j    <- x[dataPair[1]]
    g_k    <- x[dataPair[2]]
    bar_gj <- meanExpression[dataPair[1]]
    bar_gk <- meanExpression[dataPair[2]]
    return((g_j-bar_gj)*(g_k-bar_gk))
  }
  
  # calculate dataPair covariance for each sample (row) in fullData 
  return(as.numeric(apply(fullData, MARGIN = 1, FUN = coVar_calc)))
}