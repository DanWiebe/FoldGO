#' Abstract S4 class for FuncAnnotGroups object
#'
#' @slot genegroups GeneGroups. - object of GeneGroups class
#' @slot bggenes character. - vector contains background set of genes
#' @slot resultlist list. - list with filenames as keys and annotaton data frames as values
#' @slot padjmethod character. - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'                        Benjamini-Hochberg by default
#' @slot qitborder numeric. - minimal number of genes annotated to a term (1 by default)
#' @slot wholeintname character. - name of the DEG interval (initial set of genes)
#'
#'
setClass(

  "FuncAnnotGroups",

  slots = c(
    genegroups = "GeneGroups",
    bggenes = "character",
    resultlist = "list",
    padjmethod = "character",
    qitborder = "numeric",
    wholeintname = "character"
  ),

  validity = function(object) {
    if (!object@padjmethod %in% p.adjust.methods) {
      return("Unknown multiple testing correction method.
             Choose adjustment method from p.adjust.methods")
    }
    return(TRUE)
  },

  prototype = list(
    bggenes = character(0),
    padjmethod = "BH",
    qitborder = 1
  )

)


#-------------------------Annotation Object for topGO---------------------------
#' S4 class for FuncAnnotGroupsTopGO object
#'
#' @slot namespace character. - character string specifing GO namespace ("BP", "MF" or "CC")
#' @slot genesannot numeric. - minimal number of genes annotated to a term in the annotation
#' @slot algorithm character. - from TopGO manual: character string specifing which algorithm to use. The algorithms are shown by the topGO whichAlgorithms() function.
#' @slot statistic character. - from TopGO manual: character string specifing which test to use. The statistical tests are shown by the topGO whichTests() function.
#' @slot annot function. - from TopGO manual: These functions are used to compile a list of GO terms such that each element in the list is a character vector containing all the gene identifiers that are mapped to the respective GO term.
#' @slot gene2GO list. - from TopGO manual: named list of character vectors.  The list names are genes identifiers.  For each gene the character vector contains the GO identifiers it maps to.  Only the most specific annotations are required.
#' @slot mapping character. - from TopGO manual: character string specifieng the name of the Bioconductor package containing the gene mappings for a specific organism. For example: mapping = "org.Hs.eg.db".
#' @slot ID character. - from TopGO manual: character string specifing the gene identifier to use. Currently only the following identifiers can be used: c("entrez", "genbank", "alias", "ensembl", "symbol",     "genename", "unigene")
#'
#'
setClass(

  "FuncAnnotGroupsTopGO",

  slots = c(
    namespace = "character",
    genesannot = "numeric",
    algorithm = "character",
    statistic = "character",
    annot = "function",
    gene2GO = "list",
    mapping = "character",
    ID = "character"
  ),

  prototype = list(
    genesannot = 1,
    algorithm = "classic",
    statistic = "fisher",
    annot = topGO::annFUN,
    gene2GO = list(0),
    mapping = "custom",
    ID = character(0)
  ),

  validity = function(object){
    if (!object@namespace %in% c("BP", "MF", "CC")) {
      return("Unknown namespace. Choose one of the following: MF, CC, BP")
    }
    return(TRUE)
  },

  contains = "FuncAnnotGroups"

)

#' Constructor for FuncAnnotGroupsTopGO S4 class
#'
#' @param genegroups - object of GeneGroups class
#' @param namespace - character string specifing GO namespace ("BP", "MF" or "CC")
#' @param ... - Other parameters:
#' \itemize{
#' \item genesannot - minimal number of genes annotated to a term in the annotation
#' \item algorithm - from TopGO manual: character string specifing which algorithm to use. The algorithms are shown by the topGO whichAlgorithms() function.
#' \item statistic - from TopGO manual: character string specifing which test to use. The statistical tests are shown by the topGO whichTests() function.
#' \item annot - from TopGO manual: These functions are used to compile a list of GO terms such that each element in the list is a character vector containing all the gene identifiers that are mapped to the respective GO term.
#' \item gene2GO - from TopGO manual: named list of character vectors.  The list names are genes identifiers.  For each gene the character vector contains the GO identifiers it maps to.  Only the most specific annotations are required.
#' \item mapping - from TopGO manual: character string specifieng the name of the Bioconductor package containing the gene mappings for a specific organism. For example: mapping = "org.Hs.eg.db".
#' \item ID - from TopGO manual: character string specifing the gene identifier to use. Currently only the following identifiers can be used: c("entrez", "genbank", "alias", "ensembl", "symbol",     "genename", "unigene").
#' }
#'
#' @description Constructor function that creates object of FuncAnnotGroupsTopGO class.
#'              It takes GeneGroups object and GO namespace ("BP", "MF" or "CC") as a minimal set of input parameters.
#'              For more details see Arguments section.
#'
#' @export
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair", package = "FsgorS4package", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path)
#' gaf_list <- convertToList(gaf)
#' annotobj <- FuncAnnotGroupsTopGO(up_groups,"BP", gene2GO = gaf_list, annot = topGO::annFUN.GO2genes, bggenes = bggenes, padjmethod = "BH", qitborder = 10, genesannot = 1)
FuncAnnotGroupsTopGO <- function(genegroups, namespace, ...) {

  if (!requireNamespace("TopGO", quietly = TRUE)) {
    stop("TopGO package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  obj <- new(
    "FuncAnnotGroupsTopGO",
    genegroups = genegroups,
    namespace = namespace,
    ...
  )

  obj@wholeintname <- getWholeIntName(genegroups)

  obj <- runFuncAnnotTest(obj)
  return(obj)
}
