#######################################################
# FakeaAccessor method which returns the same list
# which it takes as input in the case when annotation
# already presented as s GO iD - gene ID annotation list
#######################################################
#' @rdname getAnnotation
#' @export
setMethod(f = "getAnnotation",
          signature = "list",
          definition = function(object){
            return(object)
          }
)

#######################################################
# Accessor method for MgsaSets sets slot which contains
# GO iD - gene ID annotation list
#######################################################
#' @rdname getAnnotation
#' @export
setMethod(f = "getAnnotation",
          signature = "MgsaSets",
          definition = function(object){
            return(object@sets[-1])
          }
)

#######################################################
# If no custom annotation passed
#######################################################
#' @rdname getAnnotation
#' @export
setMethod(f = "getAnnotation",
          signature = "NULL" ,
          definition = function(object){
            "NULL"}
)
