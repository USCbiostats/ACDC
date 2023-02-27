#' ACDC
#'
#' @description Function for association of covariance to detect differential co-expression analysis of modules of correlated genes or molecular features and external data with modules to be defined by Partition.
#' 
#' @param fullData data frame or matrix with samples as rows, all features as columns
#' @param ILC information loss criterion for Partition, or the minimum intra-class correlation required for features to be condensed; 0 <= ILC <= 1; default is 0.50
#' @param externalVar data frame, matrix, or vector containing external variable data to be used for CCA, rows are samples
#' @param identifierList row vector of identifiers corresponding to columns in fullData (ex: Hugo symbols for genes); optional
#' @return Data frame, sorted by ascending BH FDR value, with columns 
#' 
#' \describe{
#' \item{moduleNum}{module number}
#' \item{features}{list of identifiers from input parameter "identifierList" for all features in the module; 0 if identifierList not included}
#' \item{CCA_corr}{list of CCA canonical correlation coefficients}
#' \item{CCA_pval}{Wilks-Lamda F-test p-value; t-test p-value if there are only 2 features in the module and a single external variable}
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
ACDC <- function(fullData, ILC = 0.50, externalVar, identifierList) {
  
  ## function to suppress output
  # use for p.asym
  hush = function(code){
    sink("/dev/null")
    tmp = code
    sink()
    return(tmp)
  }
  
  # ensure correct data types
  fullData <- as.data.frame(fullData)
  externalVar <- as.data.frame(externalVar)
  
  # partition and find modules
  part    <- partition(fullData, threshold = ILC)
  modules <- part$mapping_key[which(grepl("reduced_var_", part$mapping_key$variable)), ]
  
  df <- data.frame(moduleNum = as.numeric(gsub("reduced_var_", "", modules$variable)),
                   features = numeric(nrow(modules)),
                   CCA_corr = numeric(nrow(modules)),
                   CCA_pval = numeric(nrow(modules)))
  
  # populate dataframe for each module
  for (i in 1:nrow(modules)) {
    # map module numbers to genes
    df$features[i] <- list(identifierList[unlist(modules[i,4][[1]])])
    
    # calculate connectivity
    connectivity <- (combn(x = modules$indices[i][[1]],
                           m = 2, # pick pairs
                           FUN = coVar, fullData=fullData))
    
    # run CCA, save out correlation coefficients and wilks-lambda test
    cca_results     <- cancor(connectivity, externalVar, ycenter = F)
    df$CCA_corr[i]  <- list(cca_results$cor)
    if (ncol(connectivity) > 1 | ncol(externalVar) > 1) {
      df$CCA_pval[i]  <- hush(p.asym(rho = cca_results$cor,
                                     N = dim(connectivity)[1],
                                     p = dim(connectivity)[2],
                                     q = dim(externalVar)[2],
                                     tstat = "Wilks")$p.value[1])
    } else { ## if both connectivity and externalVar are one dimensional, use simple, one-tailed correlation test
      df$CCA_pval[i] <- pt(as.numeric(cca_results$cor)*(sqrt(nrow(df)-2/(1-as.numeric(cca_results$cor)^2))), 
                           df = nrow(df)-2, 
                           lower.tail = F)
    }
  }
  
  # FDR
  df <- df[order(df$CCA_pval), ]
  row.names(df) <- 1:nrow(df)
  df$BHFDR_qval <- p.adjust(df$CCA_pval, method="BH")
  
  return(df)
}