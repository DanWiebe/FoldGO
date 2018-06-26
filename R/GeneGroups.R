#' S4 class for Gene Groups
#'
#' @slot inputtable data.frame. - dataframe contains initial set of genes gene ID's in the first row
#'             and corresponding fold change values in the second row
#' @slot logfold logical. - TRUE if fold values are presented in log scale, otherwise is FALSE
#' @slot quannumber numeric. - number of quantiles (e.g. 2,3,4...)
#' @slot groups list. - list of gene sets for each quatile and all combinations
#' @slot intnames character. - vector of intervals names
#' @slot wholeintname character. - name of the interval containing all differentially expressed genes
#' @slot regtype character. - regulation type (up or down)
#'
#'
setClass(

  "GeneGroups",

  slots = c(
    inputtable = "data.frame",
    logfold = "logical",
    quannumber = "numeric",
    groups = "list",
    intnames = "character",
    wholeintname = "character",
    regtype = "character"
  ),

  prototype = list(
    logfold = TRUE
  )

)

#' Constructor for GeneGroups S4 class
#'
#' @param inputtable - dataframe contains initial set of genes gene ID's in the first row
#'             and corresponding fold change values in the second row
#' @param quannumber - number of quantiles (e.g. 2,3,4...)
#' @param logfold - TRUE if fold values are presented in log scale, otherwise is FALSE (TRUE by default)
#'
#' @description Constructor function that creates object of GeneGroups class.
#'              It takes dataframe with genes ID's and fold values, number of quantiles
#'              and logical variable which must set to TRUE if fold values are presented in logarithmic scale,
#'              otherwise it must be set to FALSE value (TRUE by default) as parameters.
#'
#' @export
#'
#' @examples
#' GeneGroups(degenes, 6)
#' GeneGroups(degenes, 10)
#'
GeneGroups <- function(inputtable, quannumber, ...) {
  obj <- new(
    "GeneGroups",
    inputtable = inputtable,
    quannumber = quannumber,
    ...
  )
  return(divToGroups(obj))
}
