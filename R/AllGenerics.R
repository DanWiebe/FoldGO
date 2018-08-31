# Annotation object
#--------------------GETTERS-------------------
#' Accessor method for list of functional annotation result tables
#'
#' @param object - object of FuncAnnotGroups class
#'
#' @return - list with filenames as keys and annotaton data frames as values
#' @export
#'
#' @examples
#' getResultList(up_annotobj)
setGeneric(name = "getResultList",
           def = function(object){
             standardGeneric("getResultList")
           }
)

#----------------------------------------------
#' Functional annotation test
#'
#' @param object - object of FuncAnnotGroupsTopGO class
#' @return - object of FuncAnnotGroupsTopGO class
#'
setGeneric(name = "runFuncAnnotTest",
           def = function(object){
             standardGeneric("runFuncAnnotTest")
           }
)

# AnnotationReader object
#--------------------GETTERS-------------------
#' Get version of instance
#'
#' @param object - Object of S4 class
#'
#' @return - version of instance
#' @export
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma",
#'                          package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path, geneid_col = 10)
#' getVersion(gaf)
setGeneric(name = "getVersion",
           def = function(object){
             standardGeneric("getVersion")
           }
)
#' Get annotation derived from annotation file
#'
#' @param object - Object that is instance of subclass of AnnotationReader class
#'
#' @return - dataframe contains annotation
#' @export
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma",
#'                           package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path, geneid_col = 10)
#' getAnnotation(gaf)
setGeneric(name = "getAnnotation",
           def = function(object){
             standardGeneric("getAnnotation")
           }
)

#----------------------------------------------
#' Read annotation file
#'
#' @param object - object that is instance of subclass of AnnotationReader class (e.g. GAFReader class)
#' @return - instance of subclass of AnnotationReader class
#'
setGeneric(name = "read",
           def = function(object){
             standardGeneric("read")
           }
)
#' Convert annotation to list contains GO term id's as keys and Gene ID's as values
#'
#' @param object - object that is instance of subclass of AnnotationReader class (e.g. GAFReader class)
#'
#' @return - list with GO term id's as keys and Gene ID's as values
#' @export
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma",
#'                           package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path, geneid_col = 10)
#' convertToList(gaf)
setGeneric(name = "convertToList",
           def = function(object){
             standardGeneric("convertToList")
           }
)
# FoldSpecTest object
#--------------------GETTERS-------------------
#' Get dataframe with fold-change-specific terms
#'
#' @param object - Object of FoldSpecTest class
#'
#' @return - dataframe with fold-specific GO terms data
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_annotobj)
#' getFStable(fs_up)
setGeneric(name = "getFStable",
           def = function(object){
             standardGeneric("getFStable")
           }
)
#' Get dataframe with not fold-change-specific terms
#'
#' @param object - Object of FoldSpecTest class
#'
#' @return - dataframe with not fold-specific GO terms data
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_annotobj)
#' getNFStable(fs_up)
setGeneric(name = "getNFStable",
           def = function(object){
             standardGeneric("getNFStable")
           }
)
#' Get dataframe with both fold-change-specific
#' and not fold-change-specific terms
#'
#' @param object - Object of FoldSpecTest class
#'
#' @return - table contains all GO terms and related data
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_annotobj)
#' getResultTable(fs_up)
setGeneric(name = "getResultTable",
           def = function(object){
             standardGeneric("getResultTable")
           }
)
#----------------------------------------------
#' Find fold-change-specific terms
#'
#' @param object - Object of FoldSpecTest class
#' @param fdrstep2 - FDR threshold for 2 step of fold-specificty recognition procedure
#'
#' @return object of FoldSpecTest class
#'
setGeneric(name = "findFSterms",
           def = function(object, fdrstep2 = NULL){
             standardGeneric("findFSterms")
           }
)

#' Run enrichment analysis
#'
#' @param object - Object of FoldSpecTest class
#'
#' @return object of FoldSpecTest class
#'
setGeneric(name = "calcFSsignificance",
           def = function(object){
             standardGeneric("calcFSsignificance")
           }
)
# GeneGroups generics
#--------------------GETTERS-------------------
#' Get list of gene sets for each quatile and all combinations
#'
#' @param object - Object of GeneGroups class
#'
#' @return - list of gene sets for each quatile and all combinations
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getGroups(gene_groups)
setGeneric(name = "getGroups",
           def = function(object){
             standardGeneric("getGroups")
           }
)
#' Get number of quantiles
#'
#' @param object - Object of GeneGroups class
#'
#' @return - number of quantiles
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getQuanNumber(gene_groups)
setGeneric(name = "getQuanNumber",
           def = function(object){
             standardGeneric("getQuanNumber")
           }
)
#' Get vector of intervals names
#'
#' @param object - Object of GeneGroups class
#'
#' @return - vector of intervals names
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getIntNames(gene_groups)
setGeneric(name = "getIntNames",
           def = function(object){
             standardGeneric("getIntNames")
           }
)
#' Get name of the interval containing all differentially expressed genes
#'
#' @param object - Object of GeneGroups or FoldSpecTest class
#'
#' @return - name of the interval containing all differentially expressed genes
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getWholeIntName(gene_groups)
#' fs_up <- FoldSpecTest(up_annotobj)
#' getWholeIntName(fs_up)
setGeneric(name = "getWholeIntName",
           def = function(object){
             standardGeneric("getWholeIntName")
           }
)
#' Get regulation type
#'
#' @param object - Object of GeneGroups class
#'
#' @return - regulation type (up or down)
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getRegType(gene_groups)
setGeneric(name = "getRegType",
           def = function(object){
             standardGeneric("getRegType")
           }
)
#----------------------------------------------
#' Divide initial set of genes in to quantiles
#'
#' @param object - Object of GeneGroups class
#' @return object of GeneGroups class
#'
setGeneric(name = "divToGroups",
           def = function(object){
             standardGeneric("divToGroups")
           }
)
