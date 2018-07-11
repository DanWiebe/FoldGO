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
#'
setClass(

  "GAFReader",

  slots = c(
    version = "character",
    info = "character"
  ),

  contains = "AnnotationReader"

)

#' Constructor for GAFReader S4 class
#'
#' @param file - full path to annotation file
#'
#' @description Constructor function that creates object of GAFReader class.
#'              As a parameter it takes full path to file of GAF format.
#'
#' @return - object of GAFReader class
#' @export
#'
#' @examples
#' gaf_path <- system.file("extdata", "gene_association.tair.lzma", package = "FoldGO", mustWork = TRUE)
#' gaf <- GAFReader(file = gaf_path)
GAFReader <- function(file) {
  obj <- new(
    "GAFReader",
    file = file
  )
  return(read(obj))
}
