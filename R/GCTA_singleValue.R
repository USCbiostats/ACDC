#' GCTA_singleValue
#'
#' @description Function to return the heritability of an external phenotype for a single dataset
#' 
#' @param fileLoc absolute file path to bed, bim, and fam files, including prefix
#' @param externalVar vector of length n of external variable values with no ID column; must be in the same sample order as bed, bim, fam files
#' @param gctaPath absolute path to GCTA software
#' @param remlAlg algorithm to run REML iterations in GCTA; 0 = average information (AI), 1 = Fisher-scoring, 2 = EM; default is 0 (AI)
#' @param maxRemlIt the maximum number of REML iterations; default is 100
#' @param numCovars n x c_n matrix of numerical covariates to adjust heritability model for; must be in same person order as fam file; default is NULL
#' @param catCovars n x c_c matrix of categorical covariates to adjust heritability model for; must be in same person order as fam file; default is NULL
#' @return Row of GREML output containing heritability point estimate of external data and standard error
#' 
#' @details Genome-wide Complex Trait Analysis (GCTA) is a suite of C++ functions. In order to use the GCTA functions, the user must specify the absolute path to the GCTA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/gcta/#Download).
#' 
#' Here, we use GCTA's Genomics REstricted Maximum Likelihood (GREML) method to estimate the heritability of an external phenotype.
#' 
#' @examples 
#' 
#' externalVar <- c()
#' 
#' # run function; input data before running
#' \dontrun{OSCA_singleValue(fileLoc = "pathHere", 
#'                   externalVar = externalVar,
#'                   gctaPath = "pathHere")}
#' 
#' @references 
#' 
#' Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011 Jan 7;88(1):76-82. doi: 10.1016/j.ajhg.2010.11.011. Epub 2010 Dec 17. PMID: 21167468; PMCID: PMC3014363.
#' 
#' @seealso GCTA software - \url{https://yanglab.westlake.edu.cn/software/gcta/}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
GCTA_singleValue <- function(fileLoc, 
                             externalVar, 
                             gctaPath,
                             remlAlg = 0,
                             maxRemlIt = 100,
                             numCovars = NULL, 
                             catCovars = NULL) {
  
  # check parameters
  if(!(remlAlg %in% c(0,1,2))) stop("remlAlg must be 0, 1, or 2.")
  if(!is.numeric(maxRemlIt)) stop("maxRemlIt must be numeric.")
  if(maxRemlIt < 0) stop("maxRemlIt must be positive.")
  
  ## 1. save out external data in temporary file
  # read in fam file
  fam.file <- genio::read_fam(paste0(fileLoc, ".fam"), verbose = FALSE)
  
  # write out external data
  evar_file <- data.frame(IID = fam.file$fam,
                          IID_copy = fam.file$id,
                          evar = externalVar)
  phf <- tempfile(pattern = "GCTA", fileext = ".txt")
  write.table(evar_file, phf, row.names = F, quote = F)
  
  ## 2. save out covariate data in temporary files
  # numerical covariates
  if(!is.null(numCovars)) {
    numCovars <- cbind(fam.file$fam, fam.file$id, numCovars)
    ncf       <- tempfile(pattern = "GCTA", fileext = ".txt")
    write.table(numCovars, ncf, row.names = F, col.names = F, quote = F)
  }
  
  # categorical covariates
  if(!is.null(catCovars)) {
    catCovars <- cbind(fam.file$fam, fam.file$id, catCovars)
    ccf       <- tempfile(pattern = "GCTA", fileext = ".txt")
    write.table(catCovars, ccf, row.names = F, col.names = F, quote = F)
  }
  
  ## 3. run GREML
  # create GRM
  GRM <- tempfile(pattern = "GCTA", fileext = ".txt")
  invisible(system(paste0(gctaPath, "/gcta-1.94.1 --bfile ", fileLoc, 
                          " --autosome --maf 0.1 --make-grm --out ", GRM),
                   intern=T))
  
  # estimate heritability
  herit <- tempfile(pattern = "GCTA", fileext = ".txt")
  if(!is.null(catCovars) & !is.null(numCovars)) {
    # include both numeric and categorical covariates
    invisible(system(paste0(gctaPath, "/gcta-1.94.1 --grm ", GRM, " --pheno ", phf, 
                            " --covar ", ccf, " --qcovar ", ncf, " --reml --reml-alg ", remlAlg, 
                            " --reml-maxit ", maxRemlIt, " --out ", herit),
                     intern=T))
  } else if(!is.null(catCovars)) {
    # include only categorical covariates
    invisible(system(paste0(gctaPath, "/gcta-1.94.1 --grm ", GRM, " --pheno ", phf, 
                            " --covar ", ccf, " --reml --reml-alg ", remlAlg, 
                            " --reml-maxit ", maxRemlIt, " --out ", herit),
                     intern=T))
  } else if(!is.null(numCovars)) {
    # include only numeric covariates
    invisible(system(paste0(gctaPath, "/gcta-1.94.1 --grm ", GRM, " --pheno ", phf, 
                            " --qcovar ", ncf, " --reml --reml-alg ", remlAlg, 
                            " --reml-maxit ", maxRemlIt, " --out ", herit),
                     intern=T))
  } else {
    # no covariates
    invisible(system(paste0(gctaPath, "/gcta-1.94.1 --grm ", GRM, " --pheno ", phf, 
                            " --reml --reml-alg ", remlAlg, " --reml-maxit ", maxRemlIt, " --out ", herit),
                     intern=T))
  }
  
  ## 4. save out GREML results
  data <- as.data.frame(read.csv(paste0(herit,".hsq"), sep="\t"))
  
  ## 5. delete temporary files and return results
  tmp_dir <- tempdir()
  file.remove(list.files(tmp_dir, full.names = T, pattern = "GCTA"))
  
  # to return
  data[data$Source == "V(G)/Vp", ]
}