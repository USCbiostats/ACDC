#' ACDChighdim
#'
#' @description ACDC detects differential co-expression between a set of genes, such as a module of co-expressed genes, and a set of external features (exposures or responses) by using canonical correlation analysis (CCA) on the external features and module co-expression values. A high-dimensional module is supplied by the user.
#' 
#' @param moduleIdentifier the module identifier given by Partition or other dimension reduction/clustering algorithm; default is 1
#' @param moduleCols list containing indices of column locations in fullData that specify features in the module
#' @param fullData data frame or matrix with samples as rows, all features as columns; each entry should be numeric gene expression or other molecular data values
#' @param externalVar data frame, matrix, or vector containing external variable data to be used for CCA, rows are samples; all elements must be numeric
#' @param identifierList optional row vector of identifiers, of the same length and order, corresponding to columns in fullData (ex: HUGO symbols for genes); default value is the column names from fullData
#' @param corrThreshold minimum correlation required between two features to be kept in the dataset; 0 \eqn{\leq} corrThreshold \eqn{\leq} 1; default value is 0.75
#' @return Data frame, designed to be row binded with output from other ACDC functions after removing the final column, with columns 
#' 
#' \describe{
#' \item{moduleNum}{module identifier}
#' \item{colNames}{list of column names from fullData of the features in the module}
#' \item{features}{list of identifiers from input parameter "identifierList" for all features in the module}
#' \item{CCA_corr}{list of CCA canonical correlation coefficients}
#' \item{CCA_pval}{Wilks-Lamda F-test p-value or t-test p-value}
#' \item{numPairsUsed}{number of feature pairs with correlation above corrThreshold}
#' }
#' 
#' @examples 
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # run function for diet and genotype
#' ACDChighdim(moduleIdentifier = 1,
#'             moduleCols = list(1:ncol(nutrimouse$lipid)),
#'             fullData = nutrimouse$lipid,
#'             externalVar = data.frame(diet=as.numeric(nutrimouse$diet), 
#'                                      genotype=as.numeric(nutrimouse$genotype)))
#' 
#' @details If the number of co-expression features in a particular module is larger than the number of samples, CCA will return correlation coefficients of 1, and p-values and BH FDR q-values will not be calculated. This function accepts one of these high dimension modules and reduces the dimensionality by calculating the pairwise correlations for all features and only keeping feature pairs with |correlation| > the user defined corrThreshold with a maximum number of features pairs of \eqn{\lfloor\frac{N}{2}\rfloor}. We posit that these highly correlated pairs are the skeleton structure of the full module and therefore an appropriate approximation. Once this structure is identified, co-expression values are calculated and CCA is performed as in ACDC.
#' 
#' For more information about how the co-expression features are calculated, see the coVar documentation.
#' 
#' Following CCA, which determines linear combinations of the co-expression and external feature vectors that maximize the cross-covariance matrix for each module, a Wilks-Lambda test is performed to determine if the correlation between these linear combinations is significant. If they are significant, that implies there is differential co-expression. If there is only one co-expression value for a module (ie two features in the module) and a single external variable, CCA reduces to a simple correlation test, and the t-distribution is used to test for significant correlation (Widmann, 2005). 
#'
#' @references 
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal statistical society: series B (Methodological)* **57** (1995) 289–300.
#' 
#' Martin P, et al. Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed through a nutrigenomic study. *Hepatology*, in press, 2007.
#' 
#' Queen K, Nguyen MN, Gilliland F, Chun S, Raby BA, Millstein J. ACDC: a general approach for detecting phenotype or exposure associated co-expression. (in press). *Frontiers in Medicine* (2023).
#' 
#' Widmann M. One-Dimensional CCA and SVD, and Their Relationship to Regression Maps. *Journal of Climate* **18** (2005) 2785–2792. doi:10.1175/jcli3424.1.
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import CCA
#' @import CCP
#' @import utils
#' @import stats
#' @import tidyr
ACDChighdim <- function(moduleIdentifier = 1, moduleCols, fullData, externalVar, identifierList=colnames(fullData), corrThreshold = 0.75) {
  
  # to remove "no visible binding" note
  moduleNum <- CCA_pval <- numPairsUsed <- NULL
  
  ## function to suppress output
  # use for p.asym
  hush = function(code){
    sink(nullfile())
    tmp = code
    sink()
    return(tmp)
  }
  
  # ensure correct data types
  fullData    <- as.data.frame(fullData)
  externalVar <- as.data.frame(externalVar)
  if(is.null(identifierList)) identifierList <- colnames(fullData)
  
  # check correct dimensions of input
  if(nrow(fullData) != nrow(externalVar)) stop("fullData and externalVar must have the same number of rows.")
  if(ncol(fullData) != length(identifierList)) stop("identifierList must be the same length as the number of columns in fullData.")
  if(length(moduleCols) == 0) stop("No modules input.")
  if(corrThreshold < 0 | corrThreshold > 1) stop("corrThreshold must be between 0 and 1.")
  
  # results set up 
  tmp    <- c(moduleIdentifier)
  tmp[2] <- list(colnames(fullData)[moduleCols[[1]]])
  tmp[3] <- list(identifierList[unlist(moduleCols)])
  
  # calculate correlation, ignoring sign
  corrMat  <- abs(cor(fullData[,unlist(moduleCols)]))
  
  # create dataframe of pairwise correlations in descending order
  vars     <- t(combn(colnames(corrMat), 2))
  colpairs <- data.frame(vars, corrMat[vars])
  colpairs <- colpairs[order(colpairs$corrMat.vars., decreasing = T), ]
  
  # subset colpairs based on user threshold
  colpairs <- subset(colpairs, colpairs$corrMat.vars. >= corrThreshold)
  
  # if still too many pairs, subset to first n/2 pairs
  if (nrow(colpairs) > floor(nrow(fullData)/2)) {
    message(nrow(colpairs), " feature pairs selected for module ", moduleIdentifier, " using a correlation threshold of ", corrThreshold,". Using the N/2 = ", floor(nrow(fullData)/2)," most correlated pairs instead.")
    
    colpairs <- colpairs[1:floor(nrow(fullData)/2), ]
  }
  
  # stop if no data pairs
  if(nrow(colpairs) == 0) {
    # set remaining values to null
    tmp[4] <- NA
    tmp[5] <- NA
    tmp[6] <- nrow(colpairs)
    
    # return same output that will rowbind with ACDC output
    results <- as.data.frame(t(tmp))
    colnames(results) <- c("moduleNum", "colNames", "features", "CCA_corr", "CCA_pval", "numPairsUsed")
    
    # unnest columns that don't need to be lists
    results <- unnest(results, c(moduleNum, numPairsUsed, CCA_pval))
    
    # tell user if no significantly correlated pairs found
    message("No pairs detected above correlation threshold of ", corrThreshold, " for module ", moduleIdentifier, ". Choose a lower threshold.")
    
    return(results)
  }
  
  # helper function to calculate covariance
  connectivity_calc <- function(x, fd) {
    return(coVar(dataPair = c(grep(x[1], colnames(fd)),
                              grep(x[2], colnames(fd))),
                 fullData = fd))
  }
  
  # calculate connectivity for each feature pair (row) in colpairs 
  connectivity <- apply(colpairs, MARGIN = 1, FUN = connectivity_calc, fd = fullData)
  
  # run CCA and calculate p-value
  cca_results <- cancor(connectivity, externalVar, ycenter = F)
  tmp[4]      <- list(cca_results$cor)
  if (ncol(connectivity) > 1 | ncol(externalVar) > 1) {
    tmp[5] <- hush(p.asym(rho = cca_results$cor,
                          N = dim(connectivity)[1],
                          p = dim(connectivity)[2],
                          q = dim(externalVar)[2],
                          tstat = "Wilks")$p.value[1])
  } else { ## if both connectivity and externalVar are one dimensional, use simple, one-tailed correlation test
    tmp[5] <- pt(as.numeric(cca_results$cor)*(sqrt(length(moduleCols)-2/(1-as.numeric(cca_results$cor)^2))), 
                 df = length(moduleCols)-2, 
                 lower.tail = F)
  }
  
  # save out number of pairs used
  tmp[6] <- nrow(colpairs)
  
  # return same output that will rowbind with ACDC output
  results <- as.data.frame(t(tmp))
  colnames(results) <- c("moduleNum", "colNames", "features", "CCA_corr", "CCA_pval", "numPairsUsed")
  
  # unnest columns that don't need to be lists
  results <- unnest(results, c(moduleNum, CCA_pval, numPairsUsed))
  
  return(results)
}
