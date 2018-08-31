#--------------------GETTERS-------------------
#' @describeIn getGroups Get list of gene sets for each
#' quatile and all combinations
#' @export
setMethod(f = "getGroups",
          signature = "GeneGroups",
          definition = function(object){
            return(object@groups)
          }
)

#' @describeIn getQuanNumber Get number of quantiles
#' @export
setMethod(f = "getQuanNumber",
          signature = "GeneGroups",
          definition = function(object){
            return(object@quannumber)
          }
)

#' @describeIn getIntNames Get vector of intervals names
#' @export
setMethod(f = "getIntNames",
          signature = "GeneGroups",
          definition = function(object){
            return(object@intnames)
          }
)

#' @describeIn getWholeIntName Get name of the interval
#' containing all differentially expressed genes
#' @export
setMethod(f = "getWholeIntName",
          signature = "GeneGroups",
          definition = function(object){
            return(object@wholeintname)
          }
)

#' @describeIn getRegType Get regulation type
#' @export
setMethod(f = "getRegType",
          signature = "GeneGroups",
          definition = function(object){
            return(object@regtype)
          }
)
#----------------------------------------------

#' @describeIn divToGroups Divide initial set of genes in to quantiles
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
