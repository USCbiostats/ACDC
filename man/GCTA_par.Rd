% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GCTA_par.R
\name{GCTA_par}
\alias{GCTA_par}
\title{GCTA_par}
\usage{
GCTA_par(
  df,
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
  verbose = TRUE
)
}
\arguments{
\item{df}{n x p data frame or matrix of numeric -omics values with no ID column}

\item{ILCincrement}{float between zero and one determining interval between tested ILC values; default is 0.05}

\item{fileLoc}{absolute file path to bed, bim, and fam files, including prefix}

\item{gctaPath}{absolute path to GCTA software}

\item{remlAlg}{algorithm to run REML iterations in GCTA; 0 = average information (AI), 1 = Fisher-scoring, 2 = EM; default is 0 (AI)}

\item{maxRemlIt}{the maximum number of REML iterations; default is 100}

\item{numCovars}{n x c_n matrix of numerical covariates to adjust heritability model for; must be in same person order as fam file; default is NULL}

\item{catCovars}{n x c_c matrix of categorical covariates to adjust heritability model for; must be in same person order as fam file; default is NULL}

\item{summaryType}{one of "coexpression" or "covariance"; determines how to summarize each module}

\item{permute}{boolean value for whether or not to calculate values for a random permutation module summary; default is true}

\item{numNodes}{number of available compute nodes for parallelization; default is 1}

\item{verbose}{logical for whether or not to display progress updates; default is TRUE}
}
\value{
Data frame with columns
\describe{
\item{ILC}{the information loss criterion used for that iteration}
\item{InformationLost}{percent information lost due to data reduction}
\item{PercentReduction}{percent of variables condensed compared to unreduced data}
\item{AveVarianceExplained_Observed}{average heritability estimate for PC1 of observed summary data}
\item{OverallSD_Observed}{standard deviation of the heritability estimates for PC1 of observed summary data}
\item{VarianceExplained_Observed}{list of heritability estimates for PC1 of observed summary for all modules}
\item{SE_Observed}{list of standard errors of the heritability estimates for PC1 of observed summary data for all modules}
\item{AveVarianceExplained_Permuted}{average heritability for PC1 of permuted summary data}
\item{OverallSD_Permuted}{standard deviation of the heritability estimates for PC1 of permuted summary data}
\item{VarianceExplained_Permuted}{list of heritability estimates for PC1 of permuted summary data for all modules}
\item{SE_Permuted}{list of standard errors of the heritability estimates for PC1 of permuted summary data for all modules}
}
}
\description{
GCTA_par determines the average heritability of the first principal component of either the co-expression or covariance of gene expression modules for a range of increasingly reduced datasets. Dimension reduction is done with Partition, where features are only condensed into modules if the intraclass correlation between the features is at least the user-supplied information loss criterion (ILC), 0 <= ILC <= 1. An ILC of one returns the full dataset with no reduction, and an ILC of zero returns one module of all input features, reducing the dataset to the mean value. For each ILC value, with the number of ILCs tested determined by input parameter ILCincrement, the function returns the point estimate and standard error of the average heritability of the first principal component of the co-expression or covariance of the gene expression modules in the reduced dataset. If input parameter permute is true, the function also returns the same values for a random permutation of the first principle component of the appropriate matrix.
}
\details{
Genome-wide Complex Trait Analysis (GCTA) is a suite of C++ functions. In order to use the GCTA functions, the user must specify the absolute path to the GCTA software, which can be downloaded from the Yang Lab website \href{https://yanglab.westlake.edu.cn/software/gcta/#Download}{here}.

Here, we use GCTA's Genomics REstricted Maximum Likelihood (GREML) method to estimate the heritability of an external phenotype. GREML is called 2*number of modules for each ILC tested if permutations are included.

Dimension reduction is done with Partition, an agglomerative data reduction method which performs both feature condensation and extraction based on a user provided information loss criterion (ILC). Feature condensation into modules are only accepted if the intraclass correlation between the features is at least the ILC. The superPartition function is called if the gene expression dataset contains more than 4,000 features.
}
\examples{

# run function; input absolute path to OSCA software before running
\dontrun{GCTA_par(df = geneExpressionData, 
          ILCincrement = 0.25, 
          fileLoc = "pathHere",
          gctaPath = "pathHere",
          summaryType = "coexpression",
          permute = TRUE,
          numNodes = 1)}

}
\references{
Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. \emph{Bioinformatics} \strong{36} (2019) 676–681. doi:10.1093/bioinformatics/ btz661.

Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011 Jan 7;88(1):76-82. doi: 10.1016/j.ajhg.2010.11.011. Epub 2010 Dec 17. PMID: 21167468; PMCID: PMC3014363.
}
\seealso{
GCTA software - \url{https://yanglab.westlake.edu.cn/software/gcta/}
}
\author{
Katelyn Queen, \email{kjqueen@usc.edu}
}
