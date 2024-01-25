#' GCTA_parPlot
#'
#' @description GCTA_parPlot creates a graph of the output from the GCTA_par function, plotting average heritability of the first principal component of either co-expression or covariance of gene modules against information lost/percent reduction for both observed and permuted data.
#' 
#' @param df output from GCTA_par function with permutations
#' @param dataName string of name of data for graph labels; default is blank
#' @param summaryType one of "coexpression" or "covariance"; how modules were summarized for GCTA calculations 
#' @return ggplot object
#' 
#' @details Genome-wide Complex Trait Analysis (GCTA) is a suite of C++ functions. In order to use the GCTA functions, the user must specify the absolute path to the GCTA software, which can be downloaded from the Yang Lab website [here](https://yanglab.westlake.edu.cn/software/gcta/#Download).
#' 
#' In GCTA_par, we use GCTA's Genomics REstricted Maximum Likelihood (GREML) method to estimate the average heritability of the first principal component of either co-expression or covariance of gene modules. The produced plot shows these heritability estimates at varying levels of dataset reduction, calculated for observed data in blue and permuted data in red. An information loss value of 0 represents the unreduced dataset, and an information loss level of 100 represents the data being reduced to the average expression of all features.
#' 
#' @examples 
#' 
#' # run OSCA_par and save output; input absolute path to OSCA software before running
#' \dontrun{par <- GCTA_par(df = geneExpressionData, 
#'           ILCincrement = 0.25, 
#'           fileLoc = "pathHere",
#'           gctaPath = "pathHere",
#'           summaryType = "coexpression",
#'           permute = TRUE,
#'           numNodes = 1)}
#' 
#' # run function
#' \dontrun{GCTA_parPlot(df=par, dataName = "Example Data", summaryType = "coexpression")}
#' 
#' @references 
#' Millstein J, Battaglin F, Barrett M, Cao S, Zhang W, Stintzing S, et al. Partition: a surjective mapping approach for dimensionality reduction. *Bioinformatics* **36** (2019) 676–681. doi:10.1093/bioinformatics/ btz661.
#' 
#' Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011 Jan 7;88(1):76-82. doi: 10.1016/j.ajhg.2010.11.011. Epub 2010 Dec 17. PMID: 21167468; PMCID: PMC3014363.
#' 
#' @seealso GCTA software - \url{https://yanglab.westlake.edu.cn/software/gcta/#Overview}
#' 
#' @author Katelyn Queen, \email{kjqueen@@usc.edu}
#' 
#' @export
#' @import ggplot2
#' @importFrom tools toTitleCase
GCTA_parPlot <- function(df, dataName = "", summaryType) {
  
  # check that summaryType is correct
  if(!(summaryType == "coexpression" | summaryType == "covariance")) stop("summaryType must be either covariance or coexpression.")
  
  # to remove "no visible binding" note
  InformationLost <- AVEObs_upper <- AVEObs_lower <- AVEPerm_upper <- AVEPerm_lower <- AveVarianceExplained_Observed <- AveVarianceExplained_Permuted <- NULL
  
  # create confidence bounds
  i = 0
  for (i in 1:nrow(df)) {
    df$AVEObs_upper[i] <- min(1, (df$AveVarianceExplained_Observed[i] + df$OverallSD_Observed[i]))
    df$AVEObs_lower[i] <- max(0, (df$AveVarianceExplained_Observed[i] - df$OverallSD_Observed[i]))
    
    df$AVEPerm_upper[i] <- min(1, (df$AveVarianceExplained_Permuted[i] + df$OverallSD_Permuted[i]))
    df$AVEPerm_lower[i] <- max(0, (df$AveVarianceExplained_Permuted[i] - df$OverallSD_Permuted[i]))
  }
  
  # create and return graph
  ggplot(data=df) +
    geom_pointrange(aes(x = InformationLost, y = AveVarianceExplained_Observed,
                        ymin = AVEObs_lower, 
                        ymax = AVEObs_upper,
                        color = "Observed")) +
    geom_line(aes(x = InformationLost, y = AveVarianceExplained_Observed,
                  color="Observed"), size=0.75) +
    geom_pointrange(aes(x = InformationLost+0.33, y = AveVarianceExplained_Permuted,
                        ymin = AVEPerm_lower, 
                        ymax = AVEPerm_upper,
                        color = "Permuted")) +
    geom_line(aes(x = InformationLost+0.33, y = AveVarianceExplained_Permuted,
                  color="Permuted"), size=0.75) +
    ylim(c(0, 1)) + 
    xlab("Information Lost") +
    ylab("Heritability") +
    ggtitle(paste0("Heritability of Gene Module ", tools::toTitleCase(summaryType), 
                   "\n in ", dataName)) +
    scale_x_continuous(breaks = seq(0, 100, by=25),
                       sec.axis = sec_axis(~ .,
                                           labels = round(df$PercentReduction, 0),
                                           breaks = df$InformationLost,
                                           name = "Percent Reduction")) +
    scale_color_manual(name = paste0("Gene ", toTitleCase(summaryType)),
                       values = c( "Observed" = "steelblue", "Permuted" = "darkred"),
                       labels = c("Observed", "Permuted")) +
    theme_bw(base_family = "Times") +
    theme(plot.title = element_text(face="bold", hjust = 0.5))
}