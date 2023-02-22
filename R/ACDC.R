#' ACDC
#'
#' @description Function for association of covariance to detect differential co-expression analysis of modules of correlated genes and external data using Partition
#' 
#' @param fullData data frame or matrix with samples as rows, all probes as columns
#' @param ILC information loss criterion for Partition; 0 <= ILC <= 1; default is 0.50
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
#' @references 
#' \itemize{
#'    \item BH FDR: \url{https://www.jstor.org/stable/2346101}
#'    \item Partition: \url{https://pubmed.ncbi.nlm.nih.gov/31504178/}
#' }
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import partition
#' @import CCA
#' @import CCP
#' @import utils
#' @import stats
ACDC <- function(fullData, ILC = 0.50, externalVar, geneList) {
  
  ## function to suppress output
  # use for p.asym
  hush = function(code){
    sink("/dev/null")
    tmp = code
    sink()
    return(tmp)
  }
  
  # partition and find modules
  part    <- partition(fullData, threshold = ILC)
  modules <- part$mapping_key[which(grepl("reduced_var_", part$mapping_key$variable)), ]
  
  df <- data.frame(moduleNum = as.numeric(gsub("reduced_var_", "", modules$variable)),
                   genes = numeric(nrow(modules)),
                   CCA_corr = numeric(nrow(modules)),
                   CCA_pval = numeric(nrow(modules)))
  
  # populate dataframe for each module
  for (i in 1:nrow(modules)) {
    # map module numbers to genes
    df$genes[i] <- list(geneList[unlist(modules[i,4][[1]])])
    
    # calculate connectivity
    connectivity <- (combn(x = modules$indices[i][[1]],
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