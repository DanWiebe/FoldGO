# Annotation object
#--------------------GETTERS-------------------
# Accessor method for list of functional annotation result tables
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
setGeneric(name = "getVersion",
           def = function(object){
             standardGeneric("getVersion")
           }
)

# Get annotation derived from annotation file
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
setGeneric(name = "convertToList",
           def = function(object){
             standardGeneric("convertToList")
           }
)
# FoldSpecTest object
#--------------------GETTERS-------------------
# Get dataframe with fold-change-specific terms
setGeneric(name = "getFStable",
           def = function(object){
             standardGeneric("getFStable")
           }
)
# Get dataframe with not fold-change-specific terms
setGeneric(name = "getNFStable",
           def = function(object){
             standardGeneric("getNFStable")
           }
)
# Get dataframe with both fold-change-specific
# and not fold-change-specific terms
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
setGeneric(name = "getGroups",
           def = function(object){
             standardGeneric("getGroups")
           }
)
# Get number of quantiles
setGeneric(name = "getQuanNumber",
           def = function(object){
             standardGeneric("getQuanNumber")
           }
)
# Get vector of intervals names
setGeneric(name = "getIntNames",
           def = function(object){
             standardGeneric("getIntNames")
           }
)
# getWholeIntName S4 method
setGeneric(name = "getWholeIntName",
           def = function(object){
             standardGeneric("getWholeIntName")
           }
)
# Get regulation type
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
