\docType{class}
\name{gafreader_class}
\alias{gafreader_class}
\alias{GAFReader}
\alias{getVersion}
\alias{convertToList}
\alias{getVersion,GAFReader-method}
\title{S4 class for GAFReader object}
\description{
Parser for annotation presented in GAF file format (.gaf).
GAFReader function returns object which contains as a dataframe annotation
as it presented in initial file. Via GAFReader accessor method one can retrieve
annotations as list GO term id's as keys and Gene ID's as values and version of file (see Accessors section).
}
\section{Constructor}{

\code{GAFReader(file = gaf_path, geneid_col = 10)}, where:
\describe{
\item{}{\code{file} - full path to annotation file}
\item{}{\code{geneid_col} - index of column with Gene ID (2 by default)}
}
}

\section{Accessors}{

In the code examples below \code{object} is an object of GAFReader class
\describe{
\item{}{\code{getVersion(object)} - returns version of GAF file}
\item{}{\code{getAnnotation(object)} - returns annotation from GAF file in form of GO ids - Gene ids list}
}
}

\section{Methods}{

In the code examples below \code{object} is an object of GAFReader class
\describe{
\item{}{\code{getAnnotation(object)} - Convert annotation to list
contains GO term id's as keys and Gene ID's as values}
}
}

\examples{
# read .gaf file (in this example gaf file with annotation for \\emph{A.thaliana} is used)
# object returned by \\code{\\link{GAFReader}} can be used by
# \\code{\\link{FuncAnnotGroupsTopGO}} function.
gaf_path <- system.file("extdata", "gene_association.tair.lzma",
                         package = "FoldGO", mustWork = TRUE)
gaf <- GAFReader(file = gaf_path, geneid_col = 10)
# get version of file
getVersion(gaf)
# get annoitation in the list format contains GO term id's as keys and Gene ID's as values
getAnnotation(gaf)
}
