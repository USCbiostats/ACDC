#' OSCA_par
#'
#' @description OSCA_par determines the percent variance explained in an external variable (exposures or responses) for a range of increasingly reduced datasets. Dimension reduction is done with Partition, where features are only condensed into modules if the intraclass correlation between the features is at least the user-supplied information loss criterion (ILC), 0 <= ILC <= 1. An ILC of one returns the full dataset with no reduction, and an ILC of zero returns one module of all input features, reducing the dataset to the mean value. For each ILC value, with the number of ILCs tested determined by input parameter ILCincrement, the function returns the point estimate and standard error of the percent variance explained in the observed external variable by the reduced dataset. If input parameter permute is true, the function also returns the same values for a random permutation of the external variable. 
#' 
#' @param df n x p data frame or matrix of numeric -omics values with no ID column
#' @param externalVar vector of length n of external variable values with no ID column
#' @param ILCincrement float between zero and one determining interval between tested ILC values; default is 0.05
#' @param oscaPath absolute path to OSCA software
#' @param remlAlg  which algorithm to run REML iterations in GCTA; 0 = average information (AI), 1 = Fisher-scoring, 2 = EM; default is 0 (AI)
#' @param maxRemlIt the maximum number of REML iterations; default is 100
#' @param numCovars n x c_n matrix of numerical covariates to adjust heritability model for; must be in same person order as externalVar; default is NULL
#' @param catCovars n x c_c matrix of categorical covariates to adjust heritability model for; must be in same person order as externalVar; default is NULL
#' @param permute boolean value for whether or not to calculate values for a random permutation of the external variable; default is true
#' @param numNodes number of available compute nodes for parallelization; default is 1
#' @param verbose logical for whether or not to display progress updates; default is TRUE
#' 
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
#' @details OmicS-data-based Complex trait Analysis (OSCA) is a suite of C++ functions. In order to use the OSCA functions, the user must specify the absolute path to the OSCA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/osca/#Download).
#' 
#' Here, we use OSCA's Omics Restricted Maximum Likelihood (OREML) method to estimate the percent of variance in an external phenotype that can be explained by an omics profile, akin to heritability estimates in GWAS. OREML is called twice for each ILC tested if permutations are included.
#' 
#' Dimension reduction is done with Partition, an agglomerative data reduction method which performs both feature condensation and extraction based on a user provided information loss criterion (ILC). Feature condensation into modules are only accepted if the intraclass correlation between the features is at least the ILC. The superPartition function is called if the gene expression dataset contains more than 4,000 features.
#' 
#' @examples 
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # run function; input absolute path to OSCA software before running
#' \dontrun{OSCA_par(df = nutrimouse$gene, 
#'           externalVar = as.numeric(nutrimouse$diet),
#'           ILCincrement = 0.25, 
#'           oscaPath = "pathHere")}
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
#' @seealso OSCA software - \url{https://yanglab.westlake.edu.cn/software/osca/}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import data.table
#' @import partition
#' @import foreach
OSCA_par <- function(df, 
                     externalVar, 
                     ILCincrement = 0.05, 
                     oscaPath,
                     remlAlg = 0,
                     maxRemlIt = 100,
                     numCovars = NULL,
                     catCovars = NULL,
                     permute = TRUE,
                     numNodes = 1,
                     verbose = TRUE) {
  
  # check parameters
  if(nrow(df) != length(externalVar)) stop("fullData and externalVar must have the same number of rows.")
  if(0 > ILCincrement | 1 < ILCincrement) stop("ILCincrement must be between 0 and 1.")
  if(!(remlAlg %in% c(0,1,2))) stop("remlAlg must be 0, 1, or 2.")
  if(!is.numeric(maxRemlIt)) stop("maxRemlIt must be numeric.")
  if(maxRemlIt < 0) stop("maxRemlIt must be positive.")
  
  # iteration counter
  i <- 0
  
  ## Function to return weighted total information lost in reduced dataset
  # df -- mapping key from partition
  # return - information lost, single value
  total_info_lost <- function(df) {
    return((1-sum(df$information*lengths(df$indices))/sum(lengths(df$indices)))*100)
  }
  
  # setup
  ILClist             <- seq(from=0, to=1, by=ILCincrement)
  externalVar_permute <- sample(externalVar)
  
  # parallel set up
  numNodes   <- min(numNodes, length(ILClist)) # ensure no more workers than jobs
  my.cluster <- parallel::makeCluster(numNodes, outfile = "") # print messages to console
  doParallel::registerDoParallel(my.cluster)
  
  if(verbose) message("Starting analysis.")
  if(dim(df)[[2]] > 4000) message("Using superPartition due to more than 4,000 features.")
  
  # for each ILC value
  results <- foreach::foreach (i = 1:length(ILClist),
                      .combine = rbind,
                      .packages = c("partition", "data.table"),
                      .export = c("OSCA_singleValue")) %dopar% {
    # vector to store results
    tmp <- c(ILClist[i])
    
    # partition for given ILC; save out information lost and percent reduction
    if(dim(df)[[2]] > 4000) {
      prt <- partition::super_partition(df, threshold = ILClist[i])
    } else {
      prt <- partition::partition(df, threshold = ILClist[i])
    }
    tmp[2] <- round(total_info_lost(prt$mapping_key), 3) #infoLost
    tmp[3] <- round((1 - ncol(prt$reduced_data)/dim(df)[[2]])*100, 3) #percRed
    
    # calculate PVE for observed phenotype
    obs    <- OSCA_singleValue(df = prt$reduced_data,
                               externalVar = externalVar,
                               oscaPath = oscaPath,
                               remlAlg = remlAlg,
                               maxRemlIt = maxRemlIt,
                               catCovars = catCovars,
                               numCovars = numCovars)
    tmp[4] <- obs[,2] #PVE_obs
    tmp[5] <- obs[,3] #SE_obs
    
    if(permute == TRUE) {
      # calculate PVE for permuted phenotype
      per <- OSCA_singleValue(df = prt$reduced_data,
                              externalVar = externalVar_permute,
                              oscaPath = oscaPath,
                              remlAlg = remlAlg,
                              maxRemlIt = maxRemlIt,
                              catCovars = catCovars,
                              numCovars = numCovars)
      tmp[6] <- per[,2] #PVE_per
      tmp[7] <- per[,3] #SE_per
    }
    
    if(verbose) message(paste0("ILC = ", ILClist[i], " complete."))
    return(tmp)
  }
  
  # close connection
  parallel::stopCluster(cl = my.cluster)
  
  # column names for results df
  results <- as.data.frame(results)
  if(permute == TRUE) {
    colnames(results) <- c("ILC", "InformationLost", "PercentReduction", "VarianceExplained_Observed", "SE_Observed", "VarianceExplained_Permuted", "SE_Permuted")
  } else {
    colnames(results) <- c("ILC", "InformationLost", "PercentReduction", "VarianceExplained_Observed", "SE_Observed")
  }
  
  # to return
  results
}