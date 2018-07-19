#--------------------GETTERS-------------------
#' Get version of GAF file
#'
#' @param object - Object of GAFReader class
#'
#' @return - version of GAF file
#' @export
setMethod(
  f = "getVersion",
  signature = "GAFReader",
  definition = function(object) {
    return(object@version)
  }
)

#' Get annotation from GAF file
#'
#' @param object - Object that is instance of subclass of AnnotationReader class
#'
#' @return - dataframe contains annotation
#' @export
setMethod(
  f = "getAnnotation",
  signature = "AnnotationReader",
  definition = function(object) {
    return(object@annotation)
  }
)

#----------------------------------------------

#' Read GAF file
#'
#' @param object - Object of GAFReader class
#' @return - object of GAFReader class
#'
setMethod(
  f = "read",
  signature = "GAFReader",
  definition = function(object) {

    lines <- readLines(object@file)

    # separate header and data parts of file
    inds <- grep("^!.*", lines)
    header <- lines[inds]
    lines <- lines[-inds]

    object@version <- sub("!gaf-version: ", "",
                   header[grep("^!gaf-version: [0-9]\\.[0-9]$", header)])
    object@info <- header

    # split lines and put them in data frame with specified column names
    df <- matrix(character(0), ncol = 17, nrow = length(lines))
    df <- t(sapply(lines, function(x) strsplit(x, "\t")[[1]]))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    rownames(df) <- NULL
    object@annotation <- df
    return(object)
  }
)

#' Convert GAF format type annotation to list contains GO term id's as keys and Gene ID's as values
#'
#' @param object - object of GAFReader class
#'
#' @return - list with GO term id's as keys and Gene ID's as values
#' @export
setMethod(
  f = "convertToList",
  signature = "GAFReader",
  definition = function(object) {
    # extract GO id's and Gene id's annotated to them
    data <- object@annotation
    annotdf <-
      data.frame(
        "GOID" = data[, 5],
        "GeneID" = data[, object@geneid_col],
        stringsAsFactors = FALSE
      )
    # convert data frame to list with GO term id's as keys
    # and corresponding gene id's as values
    annotdf[["GOID"]] <- as.factor(annotdf[["GOID"]])
    outlist <- split(annotdf, annotdf[["GOID"]])
    outlist <- lapply(outlist, function(x)
      x[, -1])
    return(outlist)
  }
)

setMethod("show", "GAFReader",
          function(object)cat(paste0("Object of GAFReader class\n",
                                     "version = ", object@version, "\n",
                                     object@info, "\n"))
)
