#' GCTA_singleValue
#'
#' @description Function to return the heritability of an external phenotype for a single dataset
#' 
#' @param fileLoc absolute file path to bed, bim, and fam files, including prefix
#' @param externalVar vector of length n of external variable values with no ID column; must be in the same sample order as bed, bim, fam files
#' @param gctaPath absolute path to GCTA software
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
#' @import data.table 
#' @import utils
#' @import genio
#' ## Function to calculate heritability for complex trait
# @param fileLoc absolute file path to bed, bim, and fam files, including prefix
# @param externalVar vector of length n of external variable values with no ID column; must be in the same sample order as bed, bim, fam files
# @param gctaPath absolute path to GCTA software
# @return Row of GREML output containing percent variance explained in external data and standard error
GCTA_singleValue <- function(fileLoc, externalVar, gctaPath) {
  
  ## 1. save out external data in temporary file
  # read in fam file
  fam.file <- read_fam(paste0(fileLoc, ".fam"), verbose = FALSE)
  
  # write out external data
  evar_file <- data.frame(IID = fam.file$fam,
                          IID_copy = fam.file$id,
                          evar = externalVar)
  phf <- tempfile(pattern = "GCTA", fileext = ".txt")
  write.table(evar_file, phf, row.names = F, quote = F)
  
  ## 2. run GREML
  # create GRM
  GRM <- tempfile(pattern = "GCTA", fileext = ".txt")
  invisible(system(paste0(gctaPath, "/gcta-1.94.1 --bfile ", fileLoc, " --autosome --maf 0.1 --make-grm --out ", GRM),
                   intern=T))
  
  # estimate heritability
  herit <- tempfile(pattern = "GCTA", fileext = ".txt")
  invisible(system(paste0(gctaPath, "/gcta-1.94.1 --grm ", GRM, " --pheno ", phf, " --reml --out ", herit), intern=T))
  
  ## 3. save out GREML results
  data <- as.data.frame(read.csv(paste0(herit,".hsq"), sep="\t"))
  
  ## 4. delete temporary files and return results
  tmp_dir <- tempdir()
  file.remove(list.files(tmp_dir, full.names = T, pattern = "GCTA"))
  return(data[data$Source == "V(G)/Vp",])
}