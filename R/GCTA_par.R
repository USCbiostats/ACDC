#' GCTA_par
#'
#' @description GCTA_par determines the average heritability of the first principal component of either the co-expression or covariance of gene expression modules for a range of increasingly reduced datasets. Dimension reduction is done with Partition, where features are only condensed into modules if the intraclass correlation between the features is at least the user-supplied information loss criterion (ILC), 0 <= ILC <= 1. An ILC of one returns the full dataset with no reduction, and an ILC of zero returns one module of all input features, reducing the dataset to the mean value. For each ILC value, with the number of ILCs tested determined by input parameter ILCincrement, the function returns the point estimate and standard error of the average heritability of the first principal component of the co-expression or covariance of the gene expression modules in the reduced dataset. If input parameter permute is true, the function also returns the same values for a random permutation of the first principle component of the appropriate matrix. 
#' 
#' @param df n x p data frame or matrix of numeric -omics values with no ID column
#' @param ILCincrement float between zero and one determining interval between tested ILC values; default is 0.05
#' @param fileLoc absolute file path to bed, bim, and fam files, including prefix
#' @param gctaPath absolute path to GCTA software
#' @param remlAlg algorithm to run REML iterations in GCTA; 0 = average information (AI), 1 = Fisher-scoring, 2 = EM; default is 0 (AI)
#' @param maxRemlIt the maximum number of REML iterations; default is 100
#' @param numCovars n x c_n matrix of numerical covariates to adjust heritability model for; must be in same person order as fam file; default is NULL
#' @param catCovars n x c_c matrix of categorical covariates to adjust heritability model for; must be in same person order as fam file; default is NULL
#' @param summaryType one of "coexpression" or "covariance"; determines how to summarize each module
#' @param permute boolean value for whether or not to calculate values for a random permutation module summary; default is true
#' @param numNodes number of available compute nodes for parallelization; default is 1
#' @param verbose logical for whether or not to display progress updates; default is TRUE
#' 
#' @return Data frame with columns
# 
#' \describe{
#' \item{ILC}{the information loss criterion used for that iteration}
#' \item{InformationLost}{percent information lost due to data reduction}
#' \item{PercentReduction}{percent of variables condensed compared to unreduced data}
#' \item{AveVarianceExplained_Observed}{average heritability estimate for PC1 of observed summary data}
#' \item{OverallSD_Observed}{standard deviation of the heritability estimates for PC1 of observed summary data}
#' \item{VarianceExplained_Observed}{list of heritability estimates for PC1 of observed summary for all modules}
#' \item{SE_Observed}{list of standard errors of the heritability estimates for PC1 of observed summary data for all modules}
#' \item{AveVarianceExplained_Permuted}{average heritability for PC1 of permuted summary data}
#' \item{OverallSD_Permuted}{standard deviation of the heritability estimates for PC1 of permuted summary data}
#' \item{VarianceExplained_Permuted}{list of heritability estimates for PC1 of permuted summary data for all modules}
#' \item{SE_Permuted}{list of standard errors of the heritability estimates for PC1 of permuted summary data for all modules}
#' }
#' 
#' @details Genome-wide Complex Trait Analysis (GCTA) is a suite of C++ functions. In order to use the GCTA functions, the user must specify the absolute path to the GCTA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/gcta/#Download).
#' 
#' Here, we use GCTA's Genomics REstricted Maximum Likelihood (GREML) method to estimate the heritability of an external phenotype. GREML is called 2*number of modules for each ILC tested if permutations are included.
#' 
#' Dimension reduction is done with Partition, an agglomerative data reduction method which performs both feature condensation and extraction based on a user provided information loss criterion (ILC). Feature condensation into modules are only accepted if the intraclass correlation between the features is at least the ILC. The superPartition function is called if the gene expression dataset contains more than 4,000 features.
#' 
#' @examples 
#' 
#' # run function; input absolute path to OSCA software before running
#' \dontrun{GCTA_par(df = geneExpressionData, 
#'           ILCincrement = 0.25, 
#'           fileLoc = "pathHere",
#'           gctaPath = "pathHere",
#'           summaryType = "coexpression",
#'           permute = TRUE,
#'           numNodes = 1)}
#' 
#' @references 
#' Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. *Bioinformatics* **36** (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
#' 
#' Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011 Jan 7;88(1):76-82. doi: 10.1016/j.ajhg.2010.11.011. Epub 2010 Dec 17. PMID: 21167468; PMCID: PMC3014363.
#' 
#' @seealso GCTA software - \url{https://yanglab.westlake.edu.cn/software/gcta/}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import data.table
#' @import partition
#' @import genieclust
#' @import genio
#' @import foreach
GCTA_par <- function(df, 
                     ILCincrement = 0.05, 
                     fileLoc, 
                     gctaPath,
                     remlAlg = 0,
                     maxRemlIt = 100,
                     numCovars = NULL, 
                     catCovars = NULL,
                     summaryType, 
                     permute = TRUE, 
                     numNodes = 1,
                     verbose = TRUE) {
  
  # check parameters
  if(0 > ILCincrement | 1 < ILCincrement) stop("ILCincrement must be between 0 and 1.")
  if(!(summaryType == "coexpression" | summaryType == "covariance")) stop("summaryType must be either covariance or coexpression.")
  if(!(remlAlg %in% c(0,1,2))) stop("remlAlg must be 0, 1, or 2.")
  if(!is.numeric(maxRemlIt)) stop("maxRemlIt must be numeric.")
  if(maxRemlIt < 0) stop("maxRemlIt must be positive.")
  
  # iteration counter
  i <- j <- k <- m <- 0
  
  ## Function to return weighted total information lost in reduced dataset
  # df -- mapping key from partition
  # @return - information lost, single value
  total_info_lost <- function(df) {
    return((1-sum(df$information*lengths(df$indices))/sum(lengths(df$indices)))*100)
  }
  
  # setup
  ILClist <- seq(from=0, to=1, by=ILCincrement)
  
  if(verbose) message("Starting analysis.")
  if(dim(df)[[2]] > 4000) message("Using superPartition due to more than 4,000 features.")
  
  # parallel set up
  numNodes   <- min(numNodes, length(ILClist)) # ensure no more workers than jobs
  my.cluster <- parallel::makeCluster(numNodes, outfile = "") # write output from each worker to console
  doParallel::registerDoParallel(my.cluster)
  
  # for each ILC value
  results <- foreach::foreach (i = 1:length(ILClist),
                      .combine = rbind,
                      .packages = c("partition", "data.table", "genieclust", "genio"),
                      .export = c("GCTA_singleValue")) %dopar% {
                        # vector to store results
                        tmp <- c(ILClist[i])
                        
                        # partition for given ILC; save out information lost, percent reduction, and modules
                        if(dim(df)[[2]] > 4000) {
                          prt <- partition::super_partition(df, threshold = ILClist[i])
                        } else {
                          prt <- partition::partition(df, threshold = ILClist[i])
                        }
                        tmp[2]  <- round(total_info_lost(prt$mapping_key), 3) #infoLost
                        tmp[3]  <- round((1 - ncol(prt$reduced_data)/dim(df)[[2]])*100, 3) #percRed
                        modules <- prt$mapping_key[grep("reduced_var_", prt$mapping_key$variable), ]$indices
                        
                        if(verbose) message(paste0("Partitition completed for ILC = ", ILClist[i], ". ", length(modules), " modules identified."))
                        
                        # if no modules, stop
                        if(length(modules) == 0)  {
                          if(verbose) message(paste0("No modules identified when ILC = ", ILClist[i]))
                          
                          # set all values to 0
                          if(permute == TRUE) {
                            tmp[4] <- tmp[5] <- tmp[6] <- tmp[7] <- tmp[8] <- tmp[9] <- tmp[10] <- tmp[11] <- 0 
                          } else {
                            tmp[4] <- tmp[5] <- tmp[6] <- tmp[7] <- 0 
                          }
                          
                          if(verbose) message(paste0("ILC = ", ILClist[i], " complete."))
                          return(tmp)
                        }
                        
                        # placeholders
                        coex.PCs  <- NULL
                        covar.PCs <- NULL
                        
                        # get appropriate matrix based on summary type
                        ## n x num mod matrix of PC1 of summary measure
                        for (j in 1:length(modules)) {
                          
                          if(summaryType == "coexpression") {
                            # first PC of expression
                            pcs.coex <- prcomp(subset(df, select = modules[j][[1]]))
                            coex.PCs <- cbind(coex.PCs, pcs.coex$x[,1]) 
                          } else if(summaryType == "covariance") {
                            # calculate covariance matrix for module
                            covarMatrix <- (combn(x = modules[j][[1]],
                                                  m = 2, # pick pairs
                                                  FUN = coVar, fullData = df))
                            
                            # if more than 1M pairs, subset randomly
                            if(dim(covarMatrix)[[2]] > 1e6) {
                              covarMatrix <- covarMatrix[, sample(ncol(covarMatrix), size=1e6)]
                            }
                            
                            # first PC of covariance
                            pcs.covar <- prcomp(covarMatrix)
                            covar.PCs <- cbind(covar.PCs, pcs.covar$x[,1])
                          }
                        }
                        
                        if(verbose) message(paste0("PCs calculated for ILC = ", ILClist[i]))
                        
                        if(summaryType == "coexpression") {
                          # COEXPRESSION - calculate heritability for each module
                          coex.herit <- NULL
                          coex.SE    <- NULL
                          coex.herit.perm <- NULL
                          coex.SE.perm    <- NULL
                          for (k in 1:ncol(coex.PCs)) {
                            # observed
                            obs <- GCTA_singleValue(fileLoc = fileLoc,
                                                    externalVar = coex.PCs[, k],
                                                    gctaPath = gctaPath,
                                                    remlAlg = remlAlg,
                                                    maxRemlIt = maxRemlIt,
                                                    catCovars = catCovars,
                                                    numCovars = numCovars)
                            coex.herit <- append(coex.herit, obs[,2])
                            coex.SE    <- append(coex.SE,    obs[,3])
                            
                            # permuted
                            if(permute == TRUE) {
                              perm <- GCTA_singleValue(fileLoc = fileLoc,
                                                       externalVar = sample(coex.PCs[, k]),
                                                       gctaPath = gctaPath,
                                                       remlAlg = remlAlg,
                                                       maxRemlIt = maxRemlIt,
                                                       catCovars = catCovars,
                                                       numCovars = numCovars)
                              coex.herit.perm <- append(coex.herit.perm, perm[,2])
                              coex.SE.perm    <- append(coex.SE.perm,    perm[,3])
                            }
                            if (verbose) {
                              if (k %% 50 == 0) message(paste0("Heritability calculated for ", k, 
                                                               " of ", ncol(coex.PCs), " PCs."))
                            }
                          }
                          
                          tmp[4] <- mean(coex.herit) # average heritability for co-expression
                          tmp[5] <- mean(coex.SE)    # average SE for co-expression
                          tmp[6] <- list(coex.herit) # heritability estimates for all modules
                          tmp[7] <- list(coex.SE)    # SE for all modules
                          
                          if(permute == TRUE) {
                            tmp[8]  <- mean(coex.herit.perm) # average heritability for co-expression
                            tmp[9]  <- mean(coex.SE.perm)    # average SE for co-expression
                            tmp[10] <- list(coex.herit.perm) # heritability estimates for all modules
                            tmp[11] <- list(coex.SE.perm)    # SE for all modules
                          }
                          
                          if(verbose) message(paste0("Heritability of co-expression calculated for ILC = ", ILClist[i]))
                        } else if(summaryType == "covariance") {
                          # COVARIANCE - calculate heritability for each module
                          covar.herit <- NULL
                          covar.SE    <- NULL
                          covar.herit.perm <- NULL
                          covar.SE.perm    <- NULL
                          for (m in 1:ncol(covar.PCs)) {
                            # observed
                            obs <- GCTA_singleValue(fileLoc = fileLoc,
                                                    externalVar = covar.PCs[, m],
                                                    gctaPath = gctaPath,
                                                    remlAlg = remlAlg,
                                                    maxRemlIt = maxRemlIt,
                                                    catCovars = catCovars,
                                                    numCovars = numCovars)
                            covar.herit <- append(covar.herit, obs[,2])
                            covar.SE    <- append(covar.SE,    obs[,3])
                            
                            # permuted
                            if(permute == TRUE) {
                              perm <- GCTA_singleValue(fileLoc = fileLoc,
                                                       externalVar = sample(covar.PCs[, m]),
                                                       gctaPath = gctaPath,
                                                       remlAlg = remlAlg,
                                                       maxRemlIt = maxRemlIt,
                                                       catCovars = catCovars,
                                                       numCovars = numCovars)
                              covar.herit.perm <- append(covar.herit.perm, perm[,2])
                              covar.SE.perm    <- append(covar.SE.perm,    perm[,3])
                            }
                            if (verbose) {
                              if (m %% 50 == 0) message(paste0("Heritability calculated for ", m, 
                                                               " of ", ncol(covar.PCs), " PCs."))
                            }
                          }
                          
                          tmp[4] <- mean(covar.herit) # average heritability for covariance
                          tmp[5] <- sd(covar.herit)   # SD estimate for heritability estimates
                          tmp[6] <- list(covar.herit) # heritability estimates for all modules
                          tmp[7] <- list(covar.SE)    # SE for all modules
                          
                          if(permute == TRUE) {
                            tmp[8]  <- mean(covar.herit.perm) # average heritability for covariance
                            tmp[9]  <- sd(covar.herit.perm)   # SD estimate for heritability estimates
                            tmp[10] <- list(covar.herit.perm) # heritability estimates for all modules
                            tmp[11] <- list(covar.SE.perm)    # SE for all modules
                          }
                          
                          if(verbose) message(paste0("Heritability of covariance calculated for ILC = ", ILClist[i]))
                        }
                        
                        if(verbose) message(paste0("Iteration ", i, " of ", length(ILClist), " complete."))
                        return(tmp)
                      }
  
  # close connection
  parallel::stopCluster(cl = my.cluster)
  
  # column names for results df
  results <- as.data.frame(results)
  if(permute == TRUE) {
    colnames(results) <- c("ILC", "InformationLost", "PercentReduction", 
                           "AveVarianceExplained_Observed", "OverallSD_Observed", "VarianceExplained_Observed",
                           "SE_Observed", "AveVarianceExplained_Permuted", "OverallSD_Permuted",
                           "VarianceExplained_Permuted", "SE_Permuted")
    results$AveVarianceExplained_Permuted <- as.numeric(results$AveVarianceExplained_Permuted)
    results$OverallSD_Permuted            <- as.numeric(results$OverallSD_Permuted)
  } else {
    colnames(results) <- c("ILC", "InformationLost", "PercentReduction", 
                           "AveVarianceExplained_Observed", "OverallSD_Observed", 
                           "VarianceExplained_Observed", "SE_Observed")
  }
  
  # fix data types
  results                               <- data.frame(results)
  results$ILC                           <- as.numeric(results$ILC)
  results$InformationLost               <- as.numeric(results$InformationLost)
  results$PercentReduction              <- as.numeric(results$PercentReduction)
  results$AveVarianceExplained_Observed <- as.numeric(results$AveVarianceExplained_Observed)
  results$OverallSD_Observed            <- as.numeric(results$OverallSD_Observed)
  
  # to return
  results
}