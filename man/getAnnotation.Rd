\docType{methods}
\name{getAnnotation}
\alias{getAnnotation}
\alias{getAnnotation,AnnotationReader-method}
\alias{getAnnotation,list-method}
\alias{getAnnotation,MgsaSets-method}
\alias{getAnnotation,NULL-method}
\title{Get annotation derived from annotation file}
\usage{
getAnnotation(object)

\S4method{getAnnotation}{AnnotationReader}(object)

\S4method{getAnnotation}{list}(object)

\S4method{getAnnotation}{MgsaSets}(object)

\S4method{getAnnotation}{`NULL`}(object)
}
\arguments{
\item{object}{- Object of \code{mgsa} package \code{MgsaSets} class or
\code{FoldGO} \code{\link{GAFReader}} class}
}
\value{
list contains GO term id's as keys and Gene ID's as values
}
\description{
This method allows to retrieve annotation from \code{MgsaSets}
or \code{\link{GAFReader}} class object in form of list contains GO term id's
as keys and Gene ID's as values
}
\examples{
\dontrun{
gaf_path <- system.file("extdata", "gene_association.tair.lzma",
                           package = "FoldGO", mustWork = TRUE)
gaf <- GAFReader(file = gaf_path, geneid_col = 10)
getAnnotation(gaf)
}
}
