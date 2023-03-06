#' OSCA_parPlot
#'
#' @description Function to return a graph comparing percent variance explained in an external phenotype and information lost/percent reduction for both observed and permuted data
#' 
#' @param df output from OSCA_par function with permutations
#' @param externalVarName string of name of external variable for graph labels; default is blank
#' @param dataName string of name of data for graph labels; default is blank
#' @return ggplot object
#' 
#' @examples 
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # generate random phenotype
#' r.pheno <- rnorm(nrow(nutrimouse$gene))
#' 
#' ## random phenotype
#' # run OSCA_par and save output; input path to OSCA software
#' # r.par <- OSCA_par(df = nutrimouse$gene, 
#' #                 externalVar = r.pheno, 
#' #                 ILCincrement = 0.25,
#' #                 oscaPath = "pathHere", 
#' #                 numNodes = detectCores()-1, 
#' #                 permute = T)
#' 
#' # run function
#' #OSCA_parPlot(df=r.par, externalVarName = "Random Phenotype", dataName = "Nutritional Issue Genes")
#' 
#' ## observed external variables
#' # run OSCA_par and save output; input path to OSCA software
#' # par <- OSCA_par(df = nutrimouse$gene, 
#' #                 externalVar = as.numeric(nutrimouse$diet),
#' #                 ILCincrement = 0.25,
#' #                 oscaPath = "pathHere", 
#' #                 numNodes = detectCores()-1, 
#' #                 permute = T)
#' 
#' # run function
#' #OSCA_parPlot(df=par, externalVarName = "Diet", dataName = "Nutritional Issue Genes")
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
#' @import ggplot2
OSCA_parPlot <- function(df, externalVarName = "", dataName = "") {
  
  # to remove "no visible binding" note
  InformationLost <- SE_Observed <- VarianceExplained_Observed <- VarianceExplained_Permuted <- NULL
  
  # create and return graph
  return(ggplot(data=df) +
           geom_line(aes(x = InformationLost, y = VarianceExplained_Observed, color="Observed"), size=0.75) +
           geom_pointrange(aes(x = InformationLost, y = VarianceExplained_Observed,
                               ymin = VarianceExplained_Observed-SE_Observed, 
                               ymax = VarianceExplained_Observed+SE_Observed,
                               color = "Observed")) +
           geom_point(aes(x = InformationLost, y = VarianceExplained_Permuted, color="Permuted"), size=2) +
           geom_line(aes(x = InformationLost, y = VarianceExplained_Permuted, color="Permuted"), size=0.75) +
           xlab("Information Lost") +
           ylab("Percent Variance Explained") +
           ggtitle(paste0("Percent Variance Explained in ", 
                          externalVarName, "\n by ", dataName)) +
           scale_x_continuous(breaks = seq(0, 100, by=25),
                              sec.axis = sec_axis(~ .,
                                                  labels = round(df$PercentReduction,0),
                                                  breaks = df$InformationLost,
                                                  name = "Percent Reduction")) +
           scale_color_manual(name = externalVarName,
                              values = c( "Observed" = "steelblue", "Permuted" = "darkred"),
                              labels = c("Observed", "Permuted")) +
           theme_bw(base_family = "Times") +
           theme(plot.title = element_text(face="bold", hjust = 0.5)))
}