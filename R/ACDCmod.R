#' ACDCmod
#'
#' @description Function for association of covariance to detect differential co-expression analysis of modules of correlated genes and external data for user-defined module(s)
#' 
#' @param fullData data frame or matrix with samples as rows, all probes as columns
#' @param modules vector of lists where each list contains column numbers from fullData of genes included in module
#' @param externalVar data frame or matrix containing external variable data to be used for CCA, rows are samples
#' @param geneList row vector of gene names corresponding to columns in fullData
#' @return Data frame, sorted by ascending BH FDR value, with columns
#' 
#' \describe{
#' \item{moduleNum}{module number}
#' \item{genes}{list of Hugo gene symbols for all genes in the module}
#' \item{CCA_corr}{list of CCA canonical correlation coefficients}
#' \item{CCA_pval}{Wilks-Lamda F-test p-value}
#' \item{BHFDR_qval}{Benjamini-Hochberg false discovery rate q-value}
#' }
#' 
#' @references BH FDR: \url{https://www.jstor.org/stable/2346101}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import CCA
#' @import CCP
#' @import utils
#' @import stats
ACDCmod <- function(fullData, modules, externalVar, geneList) {
  
  ## function to suppress output
  # use for p.asym
  hush = function(code){
    sink("/dev/null")
    tmp = code
    sink()
    return(tmp)
  }
  
  df <- data.frame(moduleNum = c(1:length(modules)),
                   genes = numeric(length(modules)),
                   CCA_corr = numeric(length(modules)),
                   CCA_pval = numeric(length(modules)))
  
  # populate dataframe for each module
  for (i in 1:length(modules)) {
    # map module numbers to genes
    df$genes[i] <- list(geneList[unlist(modules[i][[1]])])
    
    # calculate connectivity
    connectivity <- (combn(x = modules[i][[1]],
                           m = 2, # pick pairs
                           FUN = coVar, fullData=fullData))
    
    # run CCA, save out correlation coefficients and wilks-lambda test
    cca_results     <- cancor(connectivity, externalVar, ycenter = F)
    df$CCA_corr[i]  <- list(cca_results$cor)
    df$CCA_pval[i]  <- hush(p.asym(rho = cca_results$cor,
                                   N = dim(connectivity)[1],
                                   p = dim(connectivity)[2],
                                   q = dim(externalVar)[2],
                                   tstat = "Wilks")$p.value[1])
  }
  
  # FDR
  df <- df[order(df$CCA_pval), ]
  row.names(df) <- 1:nrow(df)
  df$BHFDR_qval <- p.adjust(df$CCA_pval, method="BH")
  
  return(df)
}