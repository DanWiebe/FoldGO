#--------------------GETTERS-------------------
#' Title
#'
#' @param GAFReader - Object of GAFReader class
#'
#' @return - version of GAF file
#' @export
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma", package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path)
#' getVersion(gaf)
setMethod(
  f = "getVersion",
  signature = "GAFReader",
  definition = function(object) {
    return(object@version)
  }
)

#' Title
#'
#' @param AnnotationReader - Object that is instance of subclass of AnnotationReader class
#'
#' @return - dataframe contains annotation
#' @export
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma", package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path)
#' getAnnotation(gaf)
setMethod(
  f = "getAnnotation",
  signature = "AnnotationReader",
  definition = function(object) {
    return(object@annotation)
  }
)

#----------------------------------------------

#' Title
#'
#' @param GAFReader - Object of GAFReader class
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
    colnames(df) <- c("DB", "DB_Object_ID", "DB_Object_Symbol",
                      "Qualifier", "GO_ID", "DB:Reference",
                      "Evidence_Code", "With_(or)_From", "Aspect",
                      "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
                      "Taxon", "Date", "Assigned_By",
                      "Annotation_Extension", "Gene_Product_Form_ID")
    object@annotation <- df
    return(object)
  }
)

#' Convert GAF format type annotation to list contains GO term id's as keys and Gene ID's as values
#'
#' @param GAFReader - object of GAFReader class
#'
#' @return - list with GO term id's as keys and Gene ID's as values
#' @export
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma", package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path)
#' convertToList(gaf)
setMethod(
  f = "convertToList",
  signature = "GAFReader",
  definition = function(object) {
    # extract GO id's and Gene id's annotated to them
    data <- object@annotation
    annotdf <-
      data.frame(
        "GOID" = data$GO_ID,
        "GeneID" = data$DB_Object_Name,
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
