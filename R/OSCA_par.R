#' OSCA_par
#'
#' @description Function to return a dataframe of ILCs and percent variance explained in an external phenotype for a range of ILCs between zero and one.
#' 
#' @param df n x p data frame or matrix of gene expression values with no ID column
#' @param externalVar vector of length n of external variable values with no ID column
#' @param ILCincrement float between zero and one determining interval between tested ILC values; default is 0.05
#' @param oscaPath absolute path to OSCA software
#' @param numNodes number of available compute nodes for parallelization; default is 1
#' @param permute boolean value for whether or not to include permutations; default is true
#' @return Data frame with columns
#' 
#' \describe{
#' \item{ILC}{the information loss criterion used for that iteration}
#' \item{InformationLost}{percent information lost due to data reduction}
#' \item{PercentReduction}{percent of variables condensed compared to unreduced data}
#' \item{VarianceExplained_Observed}{percent variance explained in observed external variable by the data}
#' \item{SE_Observed}{standard error of the percent variance estimate for observed external variable}
#' \item{VarianceExplained_Permuted}{percent variance explained in permuted external variable by the data; only if input parameter "permute" is true}
#' \item{SE_Permuted}{standard error of the percent variance estimate for permuted external variable; only if input parameter "permute" is true}
#' }
#' 
#' @references 
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal statistical society: series B (Methodological) 57 (1995) 289–300.
#' 
#' Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. Bioinformatics 36 (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
#' 
#' @seealso OSCA software - \url{https://yanglab.westlake.edu.cn/software/osca}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import data.table
#' @import partition
#' @import parallel
#' @import foreach
#' @import doParallel
OSCA_par <- function(df, externalVar, ILCincrement = 0.05, oscaPath, numNodes = 1, permute = T) {
  
  # iteration counter
  i = 0
  
  ## Function to return weighted total information lost in reduced dataset
  # df -- mapping key from partition
  # return - information lost, single value
  total.info.lost <- function(df) {
    return((1-sum(df$information*lengths(df$indices))/sum(lengths(df$indices)))*100)
  }
  
  # setup
  ILClist             <- seq(from=0, to=1, by=ILCincrement)
  externalVar_permute <- sample(externalVar)
  
  # parallel set up
  my.cluster <- parallel::makeCluster(
    numNodes, 
    type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  
  # PVE for each value with permutations or PVE for each value without permutations
  if (permute == T) {
    # for each ILC value
    results <- foreach (i = 1:length(ILClist),
                        .combine = rbind,
                        .packages = c("partition", "data.table"),
                        .export = c("OSCA_singleValue")) %dopar% {
                          # vector to store results
                          tmp <- c(ILClist[i])
                          
                          # partition for given ILC; save out information lost and percent reduction
                          prt    <- partition(df, threshold = ILClist[i])
                          tmp[2] <- round(total.info.lost(prt$mapping_key), 3) #infoLost
                          tmp[3] <- round((1 - ncol(prt$reduced_data)/dim(df)[[2]])*100, 3) #percRed
                          
                          # calculate PVE for observed phenotype
                          obs    <- OSCA_singleValue(df = prt$reduced_data,
                                                     externalVar = externalVar,
                                                     oscaPath = oscaPath)
                          tmp[4] <- obs[,2] #PVE_obs
                          tmp[5] <- obs[,3] #SE_obs
                          
                          # calculate PVE for permuted phenotype
                          per <- OSCA_singleValue(df = prt$reduced_data,
                                                  externalVar = externalVar_permute,
                                                  oscaPath = oscaPath)
                          tmp[6] <- per[,2] #PVE_per
                          tmp[7] <- per[,3] #SE_per
                          
                          return(tmp)
                        }
    
    # close connection
    parallel::stopCluster(cl = my.cluster)
    
    # column names for results df
    results           <- as.data.frame(results)
    colnames(results) <- c("ILC", "InformationLost", "PercentReduction", "VarianceExplained_Observed", "SE_Observed", "VarianceExplained_Permuted", "SE_Permuted")
  } else {
    # for each ILC value
    results <- foreach (i = 1:length(ILClist),
                        .combine = rbind,
                        .packages = c("partition", "data.table"),
                        .export = c("OSCA_singleValue")) %dopar% {
                          # vector to store results
                          tmp <- c(ILClist[i])
                          
                          # partition for given ILC; save out information lost and percent reduction
                          prt    <- partition(df, threshold = ILClist[i])
                          tmp[2] <- round(total.info.lost(prt$mapping_key), 3) #infoLost
                          tmp[3] <- round((1 - ncol(prt$reduced_data)/dim(df)[[2]])*100, 3) #percRed
                          
                          # calculate PVE for observed phenotype
                          obs    <- OSCA_singleValue(df = prt$reduced_data,
                                                     externalVar = externalVar,
                                                     oscaPath = oscaPath)
                          tmp[4] <- obs[,2] #PVE_obs
                          tmp[5] <- obs[,3] #SE_obs
                          
                          return(tmp)
                        }
    
    # close connection
    parallel::stopCluster(cl = my.cluster)
    
    # column names for results df
    results           <- as.data.frame(results)
    colnames(results) <- c("ILC", "InformationLost", "PercentReduction", "VarianceExplained_Observed", "SE_Observed")
  }
  
  return(results)
}