#-------------------------GeneGroupsObject---------------------------

# S4 class for Gene Groups
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

# Constructor for GeneGroups S4 class
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

# Abstract S4 class for FuncAnnotGroups object
##############################PARAMS################################
# genegroups GeneGroups. - object of GeneGroups class
# bggenes character. - vector contains background set of genes
# slot resultlist list. - list with filenames as keys and annotaton
#                         data frames as values
# slot padjmethod character. - method for multiple testing correction
# qitborder numeric. - minimal number of genes annotated to a term
#                      (1 by default)
# wholeintname character. - name of the DEG interval (initial set of genes)
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
# S4 class for FuncAnnotGroupsTopGO object
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
    GO2genes = list(),
    mapping = "custom",
    ID = character(0)
  ),

  validity = function(object){
    if (!object@namespace %in% c("BP", "MF", "CC")) {
      return("Unknown namespace. Choose one of the following: MF, CC, BP")
    }

    if (object@mapping == "custom" && length(object@GO2genes) == 0) {
      return("Object or list with GOID - geneID annotatons is not provided!")
    }

    return(TRUE)
  },

  contains = "FuncAnnotGroups"

)

# Constructor for FuncAnnotGroupsTopGO S4 class
FuncAnnotGroupsTopGO <- function(genegroups, namespace, customAnnot = NULL, ...) {

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

  if (obj@mapping == "custom") {
    obj@GO2genes <- getAnnotation(customAnnot)
  }

  obj@wholeintname <- getWholeIntName(genegroups)

  obj <- runFuncAnnotTest(obj)
  return(obj)
}

#-------------------------AnnotationReaderObject---------------------

#--------Abstract Reader--------------
#  Abstract S4 class for AnnotationReader object
##############################PARAMS################################
# file character. - full path to annotation file
# annotation data.frame. - dataframe contains annotation
setClass(

  "AnnotationReader",

  slots = c(
    file = "character",
    annotation = "data.frame"
  )

)

#-------GAF format reader-------------
# S4 class for GAFReader object
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

# GAFReader class constructor
GAFReader <- function(file, ...) {
  obj <- new(
    "GAFReader",
    file = file,
    ...
  )
  return(read(obj))
}


#-------------------------FoldSpecTestObject-------------------------

# FoldSpecTest S4 class
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
