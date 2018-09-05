\docType{class}
\name{genegroups_class}
\alias{genegroups_class}
\alias{GeneGroups}
\title{S4 class for Gene Groups}
\description{
This function splits gene list into quantiles and generates all unions of neighbouring quantiles.
             It takes dataframe with genes ID's and fold values, number of quantiles
             and logical variable which must set to TRUE if fold values are presented in logarithmic scale,
             otherwise it must be set to FALSE value (TRUE by default) as parameters.
}
\section{Constructor}{

\code{GeneGroups(inputtable, quannumber, logfold)}, where:
\describe{
\item{}{\code{inputtable} - dataframe contains initial set of genes gene ID's in the first row
            and corresponding fold change values in the second row}
\item{}{\code{quannumber} - number of quantiles (e.g. 2,3,4...)}
\item{}{\code{logfold} - TRUE if fold values are presented in log scale, otherwise is FALSE}
}
}

\section{Accessors}{

In the code examples below \code{object} is an object of GeneGroups class
\describe{
\item{}{\code{getGroups(object)} - returns list of gene sets for each
quatile and all combinations}
\item{}{\code{getWholeIntName(object)} - returns name of the interval
containing all differentially expressed genes}
\item{}{\code{getQuanNumber(object)} - returns number of quantiles}
\item{}{\code{getIntNames(object)} - returns vector of intervals names}
\item{}{\code{getRegType(object)} - returns regulation type}
}
}

\examples{
# split initial gene set into quantiles
gene_groups <- GeneGroups(degenes, 6)
# get list of gene sets for each quatile and all combinations
getGroups(gene_groups)
# get name of the interval containing all differentially expressed genes
getWholeIntName(gene_groups)
# get number of quantiles
getQuanNumber(gene_groups)
# get vector of intervals names
getIntNames(gene_groups)
# get regulation type
getRegType(gene_groups)
}
