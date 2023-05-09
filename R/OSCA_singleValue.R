#' OSCA_singleValue
#'
#' @description Function to return the percent variance explained in an external phenotype for a single dataset
#' 
#' @param df n x p dataframe or matrix of numeric -omics values with no ID column
#' @param externalVar vector of length n of external variable values with no ID column
#' @param oscaPath absolute path to OSCA software
#' @return Row of OREML output containing percent variance explained in external data and standard error
#' 
#' @details OmicS-data-based Complex trait Analysis (OSCA) is a suite of C++ functions. In order to use the OSCA functions, the user must specify the absolute path to the OSCA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/osca/#Download).
#' 
#' Here, we use OSCA's Omics Restricted Maximum Likelihood (OREML) method to estimate the percent of variance in an external phenotype that can be explained by an omics profile, akin to heritability estimates in GWAS.
#' 
#' @examples 
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # run function; input absolute path to OSCA software before running
#' \dontrun{OSCA_singleValue(df = nutrimouse$gene, 
#'                   externalVar = as.numeric(nutrimouse$diet),
#'                   oscaPath = "pathHere")}
#' 
#' @references 
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal statistical society: series B (Methodological)* **57** (1995) 289–300.
#' 
#' Martin P, et al. Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed through a nutrigenomic study. *Hepatology*, in press, 2007.
#' 
#' Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. *Bioinformatics* **36** (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
#' 
#' Queen K, Nguyen MN, Gilliland F, Chun S, Raby BA, Millstein J. ACDC: a general approach for detecting phenotype or exposure associated co-expression. (in press). *Frontiers in Medicine* (2023).
#' 
#' @seealso OSCA software - \url{https://yanglab.westlake.edu.cn/software/osca/}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import data.table 
#' @import utils
OSCA_singleValue <- function(df, externalVar, oscaPath) {
  
  # check correct dimensions of input
  if(nrow(df) != length(externalVar)) stop("fullData and externalVar must have the same number of rows.")
  
  # ensure df is a dataframe
  df <- as.data.frame(df)
  
  ## 1. save out data in temporary files
  # gene expression data
  df$IID <- c(1:nrow(df)) # dummy IID column for OSCA
  df     <- df[,c(ncol(df),1:(ncol(df)-1))] # make dummy IID column first column
  pf     <- tempfile(pattern = "OSCA", fileext = ".txt")
  write.table(df, pf, row.names = F)
  
  # phenotype data
  evar_file <- data.frame(IID = c(1:length(externalVar)),
                          IID_copy = c(1:length(externalVar)),
                          evar = externalVar)
  phf <- tempfile(pattern = "OSCA", fileext = ".txt")
  write.table(evar_file, phf, row.names = F)
  
  ## 2. run OSCA
  # create BOD files
  bodFile <- tempfile(pattern = "OSCA", fileext = ".txt")
  invisible(system(paste0(oscaPath, " osca --efile ", pf, " --no-fid --gene-expression --make-bod --out ", bodFile), intern=T))
  
  # create ORM files
  ormFile <- tempfile(pattern = "OSCA", fileext = ".txt")
  invisible(system(paste0(oscaPath, " osca --befile ", bodFile, " --no-fid --make-orm --out ", ormFile), intern=T))
  
  # run OREML
  OREML <- tempfile(pattern = "OSCA", fileext = ".txt")
  invisible(system(paste0(oscaPath, " osca --reml --orm ", ormFile, " --pheno ", phf, " --no-fid --out ", OREML), intern=T))
  
  ## 3. save out OREML results
  data <- as.data.frame(read.csv(paste0(OREML,".rsq"), sep="\t"))
  
  ## 4. delete temporary files and return results
  tmp_dir <- tempdir()
  file.remove(list.files(tmp_dir, full.names = T, pattern = "OSCA"))
  return(data[data$Source == "V(O)/Vp",])
}