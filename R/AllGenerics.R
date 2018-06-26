# Annotation object
#--------------------GETTERS-------------------
setGeneric(name = "getResultList",
           def = function(object){
             standardGeneric("getResultList")
           }
)

#----------------------------------------------
setGeneric(name = "runFuncAnnotTest",
           def = function(object){
             standardGeneric("runFuncAnnotTest")
           }
)

# Annotationreader object
#--------------------GETTERS-------------------
setGeneric(name = "getVersion",
           def = function(object){
             standardGeneric("getVersion")
           }
)

setGeneric(name = "getAnnotation",
           def = function(object){
             standardGeneric("getAnnotation")
           }
)

#----------------------------------------------

setGeneric(name = "read",
           def = function(object){
             standardGeneric("read")
           }
)

setGeneric(name = "convertToList",
           def = function(object){
             standardGeneric("convertToList")
           }
)
# FoldSpecTest object
#--------------------GETTERS-------------------
setGeneric(name = "getFStable",
           def = function(object){
             standardGeneric("getFStable")
           }
)

setGeneric(name = "getNFStable",
           def = function(object){
             standardGeneric("getNFStable")
           }
)

setGeneric(name = "getResultTable",
           def = function(object){
             standardGeneric("getResultTable")
           }
)
#----------------------------------------------
setGeneric(name = "findFSterms",
           def = function(object, fdrstep2 = NULL){
             standardGeneric("findFSterms")
           }
)


setGeneric(name = "calcFSsignificance",
           def = function(object){
             standardGeneric("calcFSsignificance")
           }
)
# GeneGroups generics
#--------------------GETTERS-------------------
setGeneric(name = "getGroups",
           def = function(object){
             standardGeneric("getGroups")
           }
)

setGeneric(name = "getQuanNumber",
           def = function(object){
             standardGeneric("getQuanNumber")
           }
)

setGeneric(name = "getIntNames",
           def = function(object){
             standardGeneric("getIntNames")
           }
)

setGeneric(name = "getWholeIntName",
           def = function(object){
             standardGeneric("getWholeIntName")
           }
)

setGeneric(name = "getRegType",
           def = function(object){
             standardGeneric("getRegType")
           }
)
#----------------------------------------------

setGeneric(name = "divToGroups",
           def = function(object){
             standardGeneric("divToGroups")
           }
)
