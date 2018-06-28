#' FoldSpecTest S4 class
#'
#' @slot fstable data.frame. - dataframe with fold-specific GO terms data
#' @slot nfstable data.frame. - dataframe with not fold-specific GO terms data
#' @slot result_table data.frame. - table contains all GO terms and related data
#' @slot annotgroups FuncAnnotGroups. - Annotaion object
#' @slot wholeintname character. - name of the interval containing all differentially expressed genes
#' @slot padjmethod character. - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'                        Benjamini-Hochberg by default
#' @slot fdrstep1 numeric. - FDR threshold for 1 step of fold-specificty recognition procedure
#' @slot fdrstep2 numeric. - FDR threshold for 2 step of fold-specificty recognition procedure
#' @slot fisher_alternative character. - indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2 by 2 case.
#'
setClass(

  "FoldSpecTest",

  slots = c(
    fstable = "data.frame",
    nfstable = "data.frame",
    result_table = "data.frame",
    annotgroups = "FuncAnnotGroups",
    wholeintname = "character",
    padjmethod = "character",
    fdrstep1 = "numeric",
    fdrstep2 = "numeric",
    fisher_alternative = "character"
  ),

  validity = function(object){
    if (!object@padjmethod %in% p.adjust.methods) {
      return("Unknown multiple testing correction method.
             Choose adjustment method from p.adjust.methods")
    }

    if (!object@fisher_alternative %in% c("greater", "less", "two.sided", "g", "l", "t")){
      return("Wrong fisher exact test alternative.
             Choose one of the following: greater, less, two.sided or shortened version: g, l, t")
    }

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("tidyr package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    return(TRUE)
  },

  prototype = list(
    padjmethod = "BH",
    fdrstep1 = 1,
    fdrstep2 = 0.05,
    fisher_alternative = "greater"
  )

)

#' Constructor for FoldSpecTest S4 class
#'
#' @param annotgroups - Annotaion object
#' @param ... - Other parameters:
#' \itemize{
#' \item fdrstep1 - FDR threshold for 1 step of fold-specificty recognition procedure (1 by default)
#' \item fdrstep2 - FDR threshold for 2 step of fold-specificty recognition procedure (0.05 by default)
#' \item padjmethod - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'                    Benjamini-Hochberg ("BH") by default
#' \item fisher_alternative - indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#'                            You can specify just the initial letter. ("greater" by default)
#' }
#'
#' @description Constructor function that creates object of FoldSpecTest class.
#'              It takes object which is instance of subclass of AnnotGroups class (e.g. FuncAnnotGroupsTopGO class)
#'              as a minimal set of input parameters. For more details see Arguments section.
#'
#' @export
#'
#' @examples
#' FoldSpecTest(up_annotobj)
#' FoldSpecTest(up_annotobj, fdrstep1 = 0.2, fdrstep2 = 0.01, padjmethod = "BY")
FoldSpecTest <- function(annotgroups, ...) {
  obj <- new(
    "FoldSpecTest",
    annotgroups = annotgroups,
    wholeintname = annotgroups@wholeintname,
    ...
  )
  obj <- calcFSsignificance(obj)
  obj <- findFSterms(obj)
  return(obj)
}
