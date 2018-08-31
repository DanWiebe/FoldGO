#' Data from rna-seq experiment on auxin treatment of Arabidopsis thaliana roots
#'
#' A dataset containing the GeneIDs and corresponding fold-change values.
#'
#' @format A data frame with 18039 rows and 4 variables:
#' \describe{
#'   \item{GeneID}{Gene identifier}
#'   \item{FC}{fold-change value}
#'   \item{pval}{p-value}
#'   \item{qval}{Benjamini-Yekutieli adjusted p-value}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97258}
"rna_seq_data"
