#' OSCA_parPlot
#'
#' @description OSCA_parPlot creates a graph of the output from the OSCA_par function, plotting percent variance explained in an external variable (exposure or response) against information lost/percent reduction for both observed and permuted data.
#' 
#' @param df output from OSCA_par function with permutations
#' @param externalVarName string of name of external variable for graph labels; default is blank
#' @param dataName string of name of data for graph labels; default is blank
#' @return ggplot object
#' 
#' @details OmicS-data-based Complex trait Analysis (OSCA) is a suite of C++ functions. In order to use the OSCA functions, the user must specify the absolute path to the OSCA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/osca/#Download).
#' 
#' In OSCA_par, we use OSCA's Omics Restricted Maximum Likelihood (OREML) method to estimate the percent of variance in an external phenotype that can be explained by an omics profile, akin to heritability estimates in GWAS. The produced plot shows the percent variance explained in an external variable at varying levels of dataset reduction, calculated for observed external variables in blue and permuted external variables in red. An information loss value of 0 represents the unreduced dataset, and an information loss level of 100 represents the data being reduced to the average expression of all features.
#' 
#' @examples 
#' #load CCA package for example dataset
#' library(CCA)
#' 
#' # load dataset
#' data("nutrimouse")
#' 
#' # run OSCA_par and save output; input absolute path to OSCA software before running
#' \dontrun{par <- OSCA_par(df = nutrimouse$gene, 
#'                 externalVar = as.numeric(nutrimouse$diet),
#'                  ILCincrement = 0.25,
#'                  oscaPath = "pathHere")}
#' 
#' # run function
#' \dontrun{OSCA_parPlot(df=par, externalVarName = "Diet", dataName = "Nutritional Issue Genes")}
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