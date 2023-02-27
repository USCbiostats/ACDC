#' OSCA_singleValue
#'
#' @description Function to return the percent variance explained in an external phenotype for a single dataset
#' 
#' @param df n x p dataframe or matrix of gene expression values with no ID column
#' @param externalVar vector of length n of external variable values with no ID column
#' @param oscaPath absolute path to OSCA software
#' @return Row of OREML output containing percent variance explained in external data and standard error
#' 
#' @references 
#' \itemize{
#'    \item BH FDR: \url{https://www.jstor.org/stable/2346101}
#'    \item Partition: \url{https://pubmed.ncbi.nlm.nih.gov/31504178/}
#' }
#' @seealso OSCA software - \url{https://yanglab.westlake.edu.cn/software/osca}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import data.table 
#' @import utils
OSCA_singleValue <- function(df, externalVar, oscaPath) {
  
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