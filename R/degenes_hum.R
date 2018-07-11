#' DEGs from rna-seq experiment on mRNA differential expression in LNCaP cells
#' expressing the wild-type androgen receptor (AR-WT) or the ligand-independent AR-V7 splice variant
#'
#' A dataset containing the GeneIDs and corresponding fold-change values.
#'
#' @format A data frame with 2079 rows and 4 variables:
#' \describe{
#'   \item{GeneID}{Gene identifier}
#'   \item{FC}{fold-change value}
#'   \item{pval}{p-value}
#'   \item{qval}{Benjamini-Yekutieli adjusted p-value}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71334}
"degenes_hum"
