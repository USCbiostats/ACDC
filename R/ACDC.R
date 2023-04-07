#' ACDC
#'
#' @description ACDC detects differential co-expression between a set of genes, such as a module of co-expressed genes, and a set of external features (exposures or responses) by using canonical correlation analysis (CCA) on the external features and module co-expression values. Modules are detected via Partition.
#' 
#' @param fullData data frame or matrix with samples as rows, all features as columns; each entry should be numeric gene expression or other molecular data values
#' @param ILC information loss criterion for Partition, or the minimum intra-class correlation required for features to be condensed; 0 <= ILC <= 1; default is 0.50
#' @param externalVar data frame, matrix, or vector containing external variable data to be used for CCA, rows are samples; all elements must be numeric
#' @param identifierList optional row vector of identifiers, of the same length and order, corresponding to columns in fullData (ex: HUGO symbols for genes); default value is the column names from fullData
#' @return Data frame, sorted by ascending BH FDR value, with columns 
#' 
#' \describe{
#' \item{moduleNum}{module identifier}
#' \item{colNames}{list of column names from fullData of the features in the module}
#' \item{features}{list of identifiers from input parameter "identifierList" for all features in the module}
#' \item{CCA_corr}{list of CCA canonical correlation coefficients}
#' \item{CCA_pval}{Wilks-Lamda F-test p-value or t-test p-value}
#' \item{BHFDR_qval}{Benjamini-Hochberg false discovery rate q-value}
#' }
#' 
#' @examples 
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # generate random phenotype
#' r.pheno <- rnorm(nrow(nutrimouse$lipid))
#' 
#' # run function for random phenotype
#' ACDC(fullData = nutrimouse$lipid,
#'      ILC = 0.50, 
#'      externalVar = r.pheno)
#' 
#' # run function for diet and genotype (non-random)
#' ACDC(fullData = nutrimouse$lipid,
#'      ILC = 0.50, 
#'      externalVar = data.frame(diet=as.numeric(nutrimouse$diet), 
#'      genotype=as.numeric(nutrimouse$genotype)))
#' 
#' @details Modules are identified by Partition, an agglomerative data reduction method which performs both feature condensation and extraction based on a user provided information loss criterion (ILC). Feature condensation into modules are only accepted if the intraclass correlation between the features is at least the ILC. For more information about how the co-expression features are calculated, see the coVar documentation.
#' 
#' Following CCA, which determines linear combinations of the co-expression and external feature vectors that maximize the cross-covariance matrix for each module, a Wilks-Lambda test is performed to determine if the correlation between these linear combinations is significant. If they are significant, that implies there is differential co-expression. If there is only one co-expression value for a module (ie two features in the module) and a single external variable, CCA reduces to a simple correlation test, and the t-distribution is used to test for significant correlation (Widmann, 2005). If the number of co-expression features in a particular module is larger than the number of samples, CCA will return correlation coefficients of 1, and p-values and BH FDR q-values will not be calculated. We are working to implement high-dimensional solutions.
#'
#' @references 
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal statistical society: series B (Methodological) 57 (1995) 289–300.
#' 
#' Martin P, et al. Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed through a nutrigenomic study. Hepatology, in press, 2007.
#' 
#' Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. Bioinformatics 36 (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
#' 
#' Widmann M. One-Dimensional CCA and SVD, and Their Relationship to Regression Maps. Journal of Climate 18 (2005) 2785–2792. doi:10.1175/jcli3424.1.
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import partition
#' @import CCA
#' @import CCP
#' @import utils
#' @import stats
ACDC <- function(fullData, ILC = 0.50, externalVar, identifierList = colnames(fullData)) {
  
  ## function to suppress output
  # use for p.asym
  hush = function(code){
    sink("/dev/null")
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
  if(0 > ILC | 1 < ILC) stop("ILC must be between 0 and 1.")
  
  # partition and find modules
  part    <- partition(fullData, threshold = ILC)
  modules <- part$mapping_key[which(grepl("reduced_var_", part$mapping_key$variable)), ]
  
  # stop if no modules identified
  if(nrow(modules) == 0) stop("No modules identified.")
  
  df <- data.frame(moduleNum = as.numeric(gsub("reduced_var_", "", modules$variable)),
                   colNames = numeric(nrow(modules)),
                   features = numeric(nrow(modules)),
                   CCA_corr = numeric(nrow(modules)),
                   CCA_pval = numeric(nrow(modules)))
  
  # populate dataframe for each module
  for (i in 1:nrow(modules)) {
    # map module numbers to column names
    df$colNames[i] <- list(colnames(fullData)[unlist(modules[i,4][[1]])])
    
    # map module numbers to identifiers
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