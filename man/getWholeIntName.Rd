\docType{methods}
\name{getWholeIntName}
\alias{getWholeIntName}
\alias{getWholeIntName,FoldSpecTest-method}
\alias{getWholeIntName,GeneGroups-method}
\title{getWholeIntName S4 method}
\arguments{
\item{object}{Object of GeneGroups or FoldSpecTest class}
}
\description{
This method returns name of the interval containing
all differentially expressed genes. It can be applied to objects of
GeneGroups and FoldSpecTest classes
}
\examples{
# GeneGroups class object example
gene_groups <- GeneGroups(degenes, 6)
getWholeIntName(gene_groups)
# FoldSpecTest class object example
fs_up <- FoldSpecTest(up_annotobj)
getWholeIntName(fs_up)
}
\seealso{
\code{\link{FoldSpecTest()}} \code{\link{GeneGroups()}}
}
