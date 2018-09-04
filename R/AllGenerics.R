# Annotation object
#--------------------GETTERS-------------------
# Accessor method for list of functional annotation result tables
#' @export
setGeneric(name = "getResultList",
           def = function(object){
             standardGeneric("getResultList")
           }
)

#----------------------------------------------
# Functional annotation test
setGeneric(name = "runFuncAnnotTest",
           def = function(object){
             standardGeneric("runFuncAnnotTest")
           }
)

# AnnotationReader object
#--------------------GETTERS-------------------
# Get version of instance
#' @export
setGeneric(name = "getVersion",
           def = function(object){
             standardGeneric("getVersion")
           }
)
# Get annotation derived from annotation file
#' @export
setGeneric(name = "getAnnotation",
           def = function(object){
             standardGeneric("getAnnotation")
           }
)

#----------------------------------------------
# Read annotation file
setGeneric(name = "read",
           def = function(object){
             standardGeneric("read")
           }
)
# Convert annotation to list contains GO term id's as keys and Gene ID's as values
#' @export
setGeneric(name = "convertToList",
           def = function(object){
             standardGeneric("convertToList")
           }
)
# FoldSpecTest object
#--------------------GETTERS-------------------
# Get dataframe with fold-change-specific terms
#' @export
setGeneric(name = "getFStable",
           def = function(object){
             standardGeneric("getFStable")
           }
)
# Get dataframe with not fold-change-specific terms
#' @export
setGeneric(name = "getNFStable",
           def = function(object){
             standardGeneric("getNFStable")
           }
)
# Get dataframe with both fold-change-specific
# and not fold-change-specific terms
#' @export
setGeneric(name = "getResultTable",
           def = function(object){
             standardGeneric("getResultTable")
           }
)
#----------------------------------------------
# Find fold-change-specific terms
setGeneric(name = "findFSterms",
           def = function(object, fdrstep2 = NULL){
             standardGeneric("findFSterms")
           }
)

# Run enrichment analysis
setGeneric(name = "calcFSsignificance",
           def = function(object){
             standardGeneric("calcFSsignificance")
           }
)
# GeneGroups generics
#--------------------GETTERS-------------------
# Get list of gene sets for each quatile and all combinations
#' @export
setGeneric(name = "getGroups",
           def = function(object){
             standardGeneric("getGroups")
           }
)
# Get number of quantiles
#' @export
setGeneric(name = "getQuanNumber",
           def = function(object){
             standardGeneric("getQuanNumber")
           }
)
# Get vector of intervals names
#' @export
setGeneric(name = "getIntNames",
           def = function(object){
             standardGeneric("getIntNames")
           }
)
#' getWholeIntName S4 method
#'
#' @description This method returns name of the interval containing
#'all differentially expressed genes. It can be applied to objects of
#'GeneGroups and FoldSpecTest classes
#'
#' @param object Object of GeneGroups or FoldSpecTest class
#'
#' @seealso \code{\link{FoldSpecTest()}} \code{\link{GeneGroups()}}
#'
#' @export
#'
#' @examples
#' # GeneGroups class object example
#' gene_groups <- GeneGroups(degenes, 6)
#' getWholeIntName(gene_groups)
#' # FoldSpecTest class object example
#' fs_up <- FoldSpecTest(up_annotobj)
#' getWholeIntName(fs_up)
setGeneric(name = "getWholeIntName",
           def = function(object){
             standardGeneric("getWholeIntName")
           }
)
# Get regulation type
#' @export
setGeneric(name = "getRegType",
           def = function(object){
             standardGeneric("getRegType")
           }
)
#----------------------------------------------
# Divide initial set of genes in to quantiles
setGeneric(name = "divToGroups",
           def = function(object){
             standardGeneric("divToGroups")
           }
)
