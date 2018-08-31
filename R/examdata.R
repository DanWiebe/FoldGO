#' Background sets of genes used in examples
#'
#' @description We used genes from two datasets in examples:
#' \enumerate{
#' \item RNA-seq experiment on auxin treatment of Arabidopsis thaliana roots (degenes)
#' \item RNA-seq experiment on mRNA differential expression in LNCaP cells expressing
#' the wild-type androgen receptor (AR-WT) or the ligand-independent AR-V7 splice variant (degenes_hum)
#' }
#' @source
#' \enumerate{
#' \item A. thaliana and auxin: \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97258}
#' \item H. sapiens LNCap AR-V7: \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71334}
#' }
#' @format A vector containing the GeneIDs with 18039 and 38238 for A. thaliana and H. sapiens correspondingly
#' @name examdata_bg
NULL

#' Precompiled objects of GeneGroups, FuncAnnotGroups
#' classes used in examples
#'
#' @description
#' \describe{
#' \item{up_groups}{object of GeneGroups class compiled from up-regulated genes
#' from rna-seq experiment on auxin treatment of Arabidopsis thaliana roots}
#' \item{down_groups}{object of GeneGroups class compiled from down-regulated genes
#' from rna-seq experiment on auxin treatment of Arabidopsis thaliana roots}
#' \item{up_annotobj}{object of FuncAnnotGroups class compiled from lists of up-regulated genes
#' from rna-seq experiment on auxin treatment of Arabidopsis thaliana roots}
#' \item{up_annotobj}{object of FuncAnnotGroups class compiled from lists of down-regulated genes
#' from rna-seq experiment on auxin treatment of Arabidopsis thaliana roots}
#' }
#'
#' @format
#' \describe{
#' \item{up_groups}{object of GeneGroups class}
#' \item{down_groups}{object of GeneGroups class}
#' \item{up_annotobj}{object of FuncAnnotGroupsTopGO class}
#' \item{up_annotobj}{object of FuncAnnotGroupsTopGO class}
#' }
#'
#' @name examdata_objs
NULL

#' Differential expressed genes used in examples
#'
#' @description We used two datasets in examples:
#' \enumerate{
#' \item RNA-seq experiment on auxin treatment of Arabidopsis thaliana roots (degenes)
#' \item RNA-seq experiment on mRNA differential expression in LNCaP cells expressing
#' the wild-type androgen receptor (AR-WT) or the ligand-independent AR-V7 splice variant (degenes_hum)
#' }
#' @source
#' \enumerate{
#' \item A. thaliana and auxin: \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97258}
#' \item H. sapiens LNCap AR-V7: \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71334}
#' }
#' @format A dataframes with 4 variables and 789 and 2079 for A. thaliana and H. sapiens correspondingly,
#' where colnames are:
#' \describe{
#'   \item{GeneID}{Gene identifier}
#'   \item{FC}{fold-change value}
#'   \item{pval}{p-value}
#'   \item{qval}{Benjamini-Yekutieli adjusted p-value}
#' }
#' @name examdata_degs
NULL

#' @rdname examdata_degs
"degenes"

#' @rdname examdata_degs
"degenes_hum"

#' @rdname examdata_bg
"bggenes"

#' @rdname examdata_bg
"bggenes_hum"

#' @rdname examdata_objs
"up_groups"

#' @rdname examdata_objs
"down_groups"

#' @rdname examdata_objs
"up_annotobj"

#' @rdname examdata_objs
"down_annotobj"
