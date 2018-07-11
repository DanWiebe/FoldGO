#--------------------GETTERS-------------------
#' Title
#'
#' @param GeneGroups - Object of GeneGroups class
#'
#' @return - list of gene sets for each quatile and all combinations
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getGroups(gene_groups)
setMethod(f = "getGroups",
          signature = "GeneGroups",
          definition = function(object){
            return(object@groups)
          }
)

#' Title
#'
#' @param GeneGroups - Object of GeneGroups class
#'
#' @return - number of quantiles
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getQuanNumber(gene_groups)
setMethod(f = "getQuanNumber",
          signature = "GeneGroups",
          definition = function(object){
            return(object@quannumber)
          }
)

#' Title
#'
#' @param GeneGroups - Object of GeneGroups class
#'
#' @return - vector of intervals names
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getIntNames(gene_groups)
setMethod(f = "getIntNames",
          signature = "GeneGroups",
          definition = function(object){
            return(object@intnames)
          }
)

#' Title
#'
#' @param GeneGroups - Object of GeneGroups class
#'
#' @return - name of the interval containing all differentially expressed genes
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getWholeIntName(gene_groups)
setMethod(f = "getWholeIntName",
          signature = "GeneGroups",
          definition = function(object){
            return(object@wholeintname)
          }
)

#' Title
#'
#' @param GeneGroups - Object of GeneGroups class
#'
#' @return - regulation type (up or down)
#' @export
#'
#' @examples
#' gene_groups <- GeneGroups(degenes, 6)
#' getRegType(gene_groups)
setMethod(f = "getRegType",
          signature = "GeneGroups",
          definition = function(object){
            return(object@regtype)
          }
)
#----------------------------------------------
#' Title
#'
#' @param GeneGroups - Object of GeneGroups class
#' @return object of GeneGroups class
#'
setMethod(f = "divToGroups",
          signature = "GeneGroups",
          definition = function(object){

            data <- object@inputtable
            n <- object@quannumber

            data <- data[order(data[, 2]), ]
            folds <- data[, 2]
            genes <- data[, 1]
            percentiles <- c(1:n) / n
            borders <- stats::quantile(folds, percentiles)
            negcheck <- sum(folds < 0)
            singleintervals <- list()
            if (negcheck == 0) {
              object@regtype <- "up"
              borders <- c(folds[1], borders)
              for (i in 1:(n - 1)) {
                singleintervals[[toString(i)]] <-
                  genes[folds >= borders[i] & folds < borders[i + 1]]
              }
              singleintervals[[toString(n)]] <-
                genes[folds >= borders[length(borders) - 1]]
            } else if (negcheck == length(folds)) {
              object@regtype <- "down"
              folds <- rev(folds)
              genes <- rev(genes)
              borders <- rev(borders)
              for (i in 1:(n - 1)) {
                singleintervals[[toString(i)]] <-
                  genes[folds <= borders[i] & folds > borders[i + 1]]
              }
              singleintervals[[toString(n)]] <-
                genes[folds <= borders[length(borders)]]
            } else {
              stop(
                "input data contains both up and down regulated genes, please analyze them separately",
                call. = FALSE
              )
            }
            crossintervals <- list()
            for (i in 1:(n - 1)) {
              finalvec <- singleintervals[[toString(i)]]
              for (j in (i + 1):n) {
                name <- paste(i, j, sep = "-")
                finalvec <- c(finalvec, singleintervals[[toString(j)]])
                crossintervals[[name]] <- unlist(finalvec)
              }
            }
            intervals <- append(singleintervals, crossintervals)

            object@groups <- intervals
            object@intnames <- names(intervals)
            object@wholeintname <- paste("1", n, sep = "-")
            return(object)
}
)

setMethod("show", "GeneGroups",
          function(object)cat(paste0("Object of GeneGroups class\n",
                                     "Number of intervals = ", object@quannumber, "\n",
                                     "Log fold values = ", object@logfold, "\n",
                                     "Regulation type = ", object@regtype))
)
