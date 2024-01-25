#' ACDCmod
#'
#' @description ACDCmod detects differential co-expression between a set of genes, such as a module of co-expressed genes, and a set of external features (exposures or responses) by using canonical correlation analysis (CCA) on the external features and module co-expression values. Modules are provided by the user.
#' 
#' @param fullData data frame or matrix with samples as rows, all probes as columns; each entry should be numeric gene expression or other molecular data values
#' @param modules vector of lists where each list contains indices of column locations in fullData that specify features in each module
#' @param externalVar data frame, matrix, or vector  containing external variable data to be used for CCA, rows are samples; all elements must be numeric
#' @param identifierList optional row vector of identifiers, of the same length and order, corresponding to columns in fullData (ex: HUGO symbols for genes); default value is the column names from fullData
#' @param numNodes number of available compute nodes for parallelization; default is 1
#' @return Tibble, sorted by ascending BH FDR value, with columns
#' 
#' \describe{
#' \item{moduleNum}{module identifier}
#' \item{colNames}{list of column names from fullData of the features in the module}
#' \item{features}{list of identifiers from input parameter "identifierList" for all features in the module}
#' \item{CCA_corr}{list of CCA canonical correlation coefficients}
#' \item{CCA_pval}{Wilks-Lamda F-test p-value; t-test p-value if there are only 2 features in the module and a single external variable}
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
#' # partition dataset and save modules
#' library(partition)
#' part <- partition(nutrimouse$lipid, threshold = 0.50)
#' mods <- part$mapping_key[which(grepl("reduced_var_", part$mapping_key$variable)), ]$mapping
#' 
#' # run function for diet and genotype
#' ACDCmod(fullData = nutrimouse$lipid,
#'         modules = mods,
#'         externalVar = data.frame(diet=as.numeric(nutrimouse$diet), 
#'                                   genotype=as.numeric(nutrimouse$genotype)))
#' 
#' @details For more information about how the co-expression features are calculated, see the coVar documentation.
#' 
#' Following CCA, which determines linear combinations of the co-expression and external feature vectors that maximize the cross-covariance matrix for each module, a Wilks-Lambda test is performed to determine if the correlation between these linear combinations is significant. If they are significant, that implies there is differential co-expression. If there is only one co-expression value for a module (ie two features in the module) and a single external variable, CCA reduces to a simple correlation test, and the t-distribution is used to test for significant correlation (Widmann, 2005). If the number of co-expression features in a particular module is larger than the number of samples, CCA will return correlation coefficients of 1, and p-values and BH FDR q-values will not be calculated. See ACDChighdim for our solution.
#'
#' @references 
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal statistical society: series B (Methodological)* **57** (1995) 289–300.
#' 
#' Martin P, et al. Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed through a nutrigenomic study. *Hepatology*, in press, 2007.
#' 
#' Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. *Bioinformatics* **36** (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
#' 
#' Queen K, Nguyen MN, Gilliland F, Chun S, Raby BA, Millstein J. ACDC: a general approach for detecting phenotype or exposure associated co-expression. *Frontiers in Medicine* (2023) 10. doi:10.3389/fmed.2023.1118824.
#' 
#' Widmann M. One-Dimensional CCA and SVD, and Their Relationship to Regression Maps. *Journal of Climate* **18** (2005) 2785–2792. doi:10.1175/jcli3424.1.
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import stats
#' @import CCP
#' @import foreach
ACDCmod <- function(fullData, modules, externalVar, identifierList=colnames(fullData), numNodes = 1) {
  
  # check correct dimensions of input
  if(nrow(fullData) != nrow(externalVar)) stop("fullData and externalVar must have the same number of rows.")
  if(ncol(fullData) != length(identifierList)) stop("identifierList must be the same length as the number of columns in fullData.")
  if(length(modules) == 0) stop("No modules input.")
  
  # to remove "no visible binding" note
  moduleNum <- CCA_pval <- NULL
  
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
  
  # iteration counter
  i = 0
  
  # parallel set up
  numNodes   <- min(numNodes, length(modules))
  my.cluster <- parallel::makeCluster(numNodes, outfile = "")
  doParallel::registerDoParallel(my.cluster)
  
  # for each module...
  results <- foreach::foreach (i = 1:length(modules),
                      .combine = rbind,
                      .packages = c("stats", "CCP"),
                      .export = c("coVar")) %dopar% {
                        # vector to store results
                        tmp <- c(i)
                        
                        # map module numbers to column names
                        ## colNames
                        tmp[2] <- list(colnames(fullData)[unlist(modules[i][[1]])])
                        
                        # map module numbers to genes
                        ## features
                        tmp[3] <- list(identifierList[unlist(modules[i][[1]])])
                        
                        # don't calculate for high dimensional modules
                        if (choose(length(modules[i][[1]]), 2) > nrow(fullData)) {
                          tmp[4] = NA ## CCA_corr
                          tmp[5] = NA ## CCA_pval
                          message("CCA values not calculated for module ", i, ". Use ACDChighdim to calculate values for this module.")
                        } else {
                          # calculate connectivity
                          connectivity <- (combn(x = modules[i][[1]],
                                                 m = 2, # pick pairs
                                                 FUN = coVar, fullData=fullData))
                          
                          # run CCA, save out correlation coefficients and wilks-lambda test
                          ## CCA_corr and CCA_pval
                          cca_results <- stats::cancor(connectivity, externalVar, ycenter = F)
                          tmp[4]      <- list(cca_results$cor)
                          if (ncol(connectivity) > 1 | ncol(externalVar) > 1) {
                            tmp[5]  <- hush(CCP::p.asym(rho = cca_results$cor,
                                                   N = dim(connectivity)[1],
                                                   p = dim(connectivity)[2],
                                                   q = dim(externalVar)[2],
                                                   tstat = "Wilks")$p.value[1])
                          } else { ## if both connectivity and externalVar are one dimensional, use simple, one-tailed correlation test
                            tmp[5] <- stats::pt(as.numeric(cca_results$cor)*(sqrt(length(modules)-2/(1-as.numeric(cca_results$cor)^2))), 
                                         df = length(modules)-2, 
                                         lower.tail = F)
                          }
                        }
                        
                        return(tmp)
                      }
  # close connection
  parallel::stopCluster(cl = my.cluster)
  
  # column names for results df
  results           <- data.frame(results)
  colnames(results) <- c("moduleNum", "colNames", "features", "CCA_corr", "CCA_pval")
  rownames(results) <- results$moduleNum
  
  # unnest columns that don't need to be lists
  results <- tidyr::unnest(results, c(moduleNum, CCA_pval))
  
  # FDR
  results            <- results[order(results$CCA_pval), ]
  results$BHFDR_qval <- p.adjust(results$CCA_pval, method="BH")
  
  # to return
  results
}
