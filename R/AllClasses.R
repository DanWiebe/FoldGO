#-------------------------GeneGroupsObject---------------------------

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
#' @param ... - Other parameters:
#' \itemize{
#' \item logfold - TRUE if fold values are presented in log scale, otherwise is FALSE
#' }
#'
#' @description Constructor function that creates object of GeneGroups class.
#'              It takes dataframe with genes ID's and fold values, number of quantiles
#'              and logical variable which must set to TRUE if fold values are presented in logarithmic scale,
#'              otherwise it must be set to FALSE value (TRUE by default) as parameters.
#' @return - object of GeneGroups class
#' @export
#' @importFrom methods new
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


#-------------------------AnnotationObject---------------------------

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
#' @slot GO2genes list. - from TopGO manual: named list of character vectors.  The list names are GO identifiers. For each GO the character vector contains the genes identifiers which are mapped to it. Only the most specific annotations are required.
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
    GO2genes = "list",
    mapping = "character",
    ID = "character"
  ),

  prototype = list(
    genesannot = 1,
    algorithm = "classic",
    statistic = "fisher",
    annot = topGO::annFUN,
    GO2genes = list(0),
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
#' \item GO2genes - from TopGO manual: named list of character vectors.  The list names are genes identifiers.  For each gene the character vector contains the GO identifiers it maps to.  Only the most specific annotations are required.
#' \item mapping - from TopGO manual: character string specifieng the name of the Bioconductor package containing the gene mappings for a specific organism. For example: mapping = "org.Hs.eg.db".
#' \item ID - from TopGO manual: character string specifing the gene identifier to use. Currently only the following identifiers can be used: c("entrez", "genbank", "alias", "ensembl", "symbol",     "genename", "unigene").
#' }
#'
#' @description Constructor function that creates object of FuncAnnotGroupsTopGO class.
#'              It takes GeneGroups object and GO namespace ("BP", "MF" or "CC") as a minimal set of input parameters.
#'              For more details see Arguments section.
#' @return - object of FuncAnnotGroupsTopGO class
#' @import topGO
#' @importFrom methods new
#' @importFrom stats p.adjust
#' @export
#'
#' @examples
#' \dontrun{
#' library(topGO)
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma",
#'                          package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path, geneid_col = 10)
#' gaf_list <- convertToList(gaf)
#' annotobj <- FuncAnnotGroupsTopGO(up_groups,"BP", GO2genes = gaf_list,
#'                                  annot = topGO::annFUN.GO2genes,
#'                                  bggenes = bggenes, padjmethod = "BH",
#'                                  qitborder = 10, genesannot = 1)
#'}
FuncAnnotGroupsTopGO <- function(genegroups, namespace, ...) {

  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("topGO package needed for this function to work. Please install it.",
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

#-------------------------AnnotationReaderObject---------------------

#--------Abstract Reader--------------
#'  Abstract S4 class for AnnotationReader object
#'
#' @slot file character. - full path to annotation file
#' @slot annotation data.frame. - dataframe contains annotation
#'
setClass(

  "AnnotationReader",

  slots = c(
    file = "character",
    annotation = "data.frame"
  )

)

#-------GAF format reader-------------
#' S4 class for GAFReader object
#'
#' @slot version character. - version of GAF file
#' @slot info character. - information from GAF file header
#' @slot geneid_col numeric. - index of column with Gene ID (2 by default)
#'
setClass(

  "GAFReader",

  slots = c(
    version = "character",
    info = "character",
    geneid_col = "numeric"
  ),

  prototype = list(
    geneid_col = 2
  ),

  contains = "AnnotationReader"

)

#' Constructor for GAFReader S4 class
#'
#' @param file - full path to annotation file
#' @param ... - Other parameters:
#' \itemize{
#' \item geneid_col - index of column with Gene ID (2 by default)
#' }
#'
#' @description Constructor function that creates object of GAFReader class.
#'              As a parameter it takes full path to file of GAF format.
#'
#' @return - object of GAFReader class
#' @export
#' @importFrom methods new
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma",
#'                          package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path, geneid_col = 10)
GAFReader <- function(file, ...) {
  obj <- new(
    "GAFReader",
    file = file,
    ...
  )
  return(read(obj))
}


#-------------------------FoldSpecTestObject-------------------------

#' FoldSpecTest S4 class
#' @name foldspectest_class
#'
#' @description FoldSpecTest object calculates test on fold-specificity and stores all resulting data needed for further analysis.
#' It takes object which is instance of subclass of AnnotGroups class (e.g. FuncAnnotGroupsTopGO class) as a minimal set of input parameters.
#' For more details see Constructor section.
#'
#' @section Constructor:
#' \code{FoldSpecTest(annotobj, fdrstep1 = 1.0, fdrstep2 = 0.05, padjmethod = "BH", fisher_alternative = "greater")}, where:
#' \describe{
#' \item{}{\code{annotobj} - object of FuncAnnotGroups class}
#' \item{}{\code{fdrstep1} - FDR threshold for 1 step of fold-specificty recognition procedure}
#' \item{}{\code{fdrstep2} - FDR threshold for 2 step of fold-specificty recognition procedure}
#' \item{}{\code{padjmethod} - method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#' Benjamini-Hochberg by default}
#' \item{}{\code{fisher_alternative} - indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' You can specify just the initial letter. Only used in the 2 by 2 case.}
#' }
#' @section Accessors:
#' In the code examples below \code{object} is an object of FoldSpecTest class
#' \describe{
#' \item{}{\code{getFStable(object)} - returns dataframe with fold-change-specific terms and related data}
#' \item{}{\code{getNFStable(object)} - returns dataframe with not fold-change-specific terms and related data}
#' \item{}{\code{getResultTable(object)} - returns dataframe with both fold-change-specific and not fold-change-specific terms}
#' \item{}{\code{getWholeIntName(object)} - returns name of largest fold-change interval (DEGs interval)}
#' }
#' @examples
#' # FoldSpecTest function requires only object of FuncAnnotGroups class as a minimal set of parameters.
#' # In the example up_annotobj is an object of FuncAnnotGroups class compiled from lists of up-regulated genes
#' # from rna-seq experiment on auxin treatment of Arabidopsis thaliana roots [FoldGO::up_annotobj].
#' FoldSpecTest(up_annotobj)
#'
#' # FoldSpecTest function with custom parameters
#' fs_up <- FoldSpecTest(up_annotobj, fdrstep1 = 0.2, fdrstep2 = 0.01, padjmethod = "BY")
#'
#' # get dataframe with fold-change-specific terms
#' getFStable(fs_up)
#'
#' # get dataframe with not fold-change-specific terms
#' getNFStable(fs_up)
#'
#' # get dataframe with both fold-change-specific and not fold-change-specific terms
#' getResultTable(fs_up)
#'
#' # get name of largest fold-change interval (DEGs interval)
#' getWholeIntName(fs_up)
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

# FoldSpecTest class constructor
#' @export
#' @importFrom methods new
#' @rdname foldspectest_class
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
