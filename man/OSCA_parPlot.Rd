% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OSCA_parPlot.R
\name{OSCA_parPlot}
\alias{OSCA_parPlot}
\title{OSCA_parPlot}
\usage{
OSCA_parPlot(df, externalVarName = "", dataName = "")
}
\arguments{
\item{df}{output from OSCA_par function with permutations}

\item{externalVarName}{string of name of external variable for graph labels; default is blank}

\item{dataName}{string of name of data for graph labels; default is blank}
}
\value{
ggplot object
}
\description{
Function to return a graph comparing percent variance explained in an external phenotype and information lost/percent reduction for both observed and permuted data
}
\examples{
#load CCA package for example dataset
library(CCA)

# load dataset
data("nutrimouse")

# generate random phenotype
r.pheno <- rnorm(nrow(nutrimouse$gene))

## random phenotype
# run OSCA_par and save output; input path to OSCA software
# r.par <- OSCA_par(df = nutrimouse$gene, 
#                 externalVar = r.pheno, 
#                 ILCincrement = 0.25,
#                 oscaPath = "pathHere", 
#                 numNodes = detectCores()-1, 
#                 permute = T)

# run function
#OSCA_parPlot(df=r.par, externalVarName = "Random Phenotype", dataName = "Nutritional Issue Genes")

## observed external variables
# run OSCA_par and save output; input path to OSCA software
# par <- OSCA_par(df = nutrimouse$gene, 
#                 externalVar = as.numeric(nutrimouse$diet),
#                 ILCincrement = 0.25,
#                 oscaPath = "pathHere", 
#                 numNodes = detectCores()-1, 
#                 permute = T)

# run function
#OSCA_parPlot(df=par, externalVarName = "Diet", dataName = "Nutritional Issue Genes")

}
\references{
Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal statistical society: series B (Methodological) 57 (1995) 289–300.

Martin P, et al. Novel aspects of PPARalpha-mediated regulation of lipid and xenobiotic metabolism revealed through a nutrigenomic study. Hepatology, in press, 2007.

Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. Bioinformatics 36 (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
}
\seealso{
OSCA software - \url{https://yanglab.westlake.edu.cn/software/osca}
}
\author{
Katelyn Queen, \email{kjqueen@usc.edu}
}
