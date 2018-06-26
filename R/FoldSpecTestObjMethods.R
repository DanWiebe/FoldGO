#--------------------GETTERS-------------------
#' Title
#'
#' @param FoldSpecTest - Object of FoldSpecTest class
#'
#' @return - dataframe with fold-specific GO terms data
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_groups)
#' getFStable(fs_up)
setMethod(f = "getFStable",
          signature = "FoldSpecTest",
          definition = function(object){
            return(object@fstable)
          }
)

#' Title
#'
#' @param FoldSpecTest - Object of FoldSpecTest class
#'
#' @return - dataframe with not fold-specific GO terms data
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_groups)
#' getNFStable(fs_up)
setMethod(f = "getNFStable",
          signature = "FoldSpecTest",
          definition = function(object){
            return(object@nfstable)
          }
)

#' Title
#'
#' @param FoldSpecTest - Object of FoldSpecTest class
#'
#' @return - table contains all GO terms and related data
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_groups)
#' getResultTable(fs_up)
setMethod(f = "getResultTable",
          signature = "FoldSpecTest",
          definition = function(object){
            return(object@result_table)
          }
)
#' Title
#'
#' @param FoldSpecTest - Object of FoldSpecTest class
#'
#' @return - name of the interval containing all differentially expressed genes
#' @export
#'
#' @examples
#' fs_up <- FoldSpecTest(up_groups)
#' getWholeIntName(fs_up)
setMethod(f = "getWholeIntName",
          signature = "FoldSpecTest",
          definition = function(object){
            return(object@wholeintname)
          }
)
#----------------------------------------------
#' Title
#'
#' @param FoldSpecTest - Object of FoldSpecTest class
#' @param fdrstep2 - FDR threshold for 2 step of fold-specificty recognition procedure
#'
setMethod(f = "findFSterms",
          signature = "FoldSpecTest",
          definition = function(object, fdrstep2 = NULL){

            if (!is.null(fdrstep2)) {
              object@fdrstep2 <- fdrstep2
            }

            result_table <- object@result_table

            fstable <- result_table[result_table$padj < object@fdrstep2, ]
            nfstable <- result_table[result_table$padj >= object@fdrstep2, ]

            if (nrow(fstable) != 0) {
               rownames(fstable) <- c(1:nrow(fstable))
            }

            if (nrow(nfstable) != 0) {
               rownames(nfstable) <- c(1:nrow(nfstable))
            }

            object@fstable <- fstable
            object@nfstable <- nfstable

            return(object)
          }
)


#' Title
#'
#' @param FoldSpecTest - Object of FoldSpecTest class
#'
setMethod(f = "calcFSsignificance",
          signature = "FoldSpecTest",
          definition = function(object){

            listoftables <- getResultList(object@annotgroups)
            fdrstep1 <- object@fdrstep1
            fdrstep2 <- object@fdrstep2
            wholeintname <- object@wholeintname
            fisher_alternative <- object@fisher_alternative
            p_adjust_method <- object@padjmethod

            # extract dataframe for whole interval annotation
            # and cut it by FDR step 2 threshold
            df <- listoftables[[wholeintname]]
            if (fdrstep1 < 1) {
              df <- subset(df, df[["padj"]] < fdrstep1)
            }
            wi_pvaldf <- df[, c("GO_id",
                                "pvalue",
                                "padj")]
            # leave only GO_id, qit and qtot columns
            # and rename qit and qtot to wqit and wqtot where w means whole interval
            df <-
              df[, c("GO_id",
                     "qit",
                     "qtot")]
            names(df)[names(df) == "qit"] <- "wqit"
            names(df)[names(df) == "qtot"] <- "wqtot"

            # create empty dataframe for merging all dataframes by GO_id
            resdf <- data.frame()

            # extract interval names form list of dataframes names
            # and delete the name of whole interval
            filenames <- names(listoftables)
            filenames <- filenames[!filenames %in% wholeintname]

            # extract dataframe for each interval except whole
            # add column with interval name
            # merge data by GO_id's from whole interval data frame
            # bind all dataframes in resdf
            for (i in 1:length(filenames)) {
              data <- listoftables[[filenames[i]]]
              data$filename <- rep(filenames[i], nrow(data))
              data <- data[, c(
                "GO_id",
                "name",
                "namespace",
                "qit",
                "qtot",
                "filename"
              )]
              resdf <-
                rbind(resdf, merge(data, df, by = "GO_id"))
            }

            # create matrix with rows contains qit, qtot, wqit, wqtot for each GO term
            cutdf <-
              resdf[, c("qit", "qtot", "wqit", "wqtot")]
            mat <- matrix(unlist(cutdf), nrow = nrow(cutdf))

            # calculate Fisher exact test for each GO term
            resdf$pvalues <-
              apply(mat, 1, function(x)
                stats::fisher.test(matrix(c(
                  x[1], x[2] - x[1], x[3] - x[1], x[4] - x[2] - x[3] + x[1]
                ), nrow = 2),
                alternative = fisher_alternative)$p.value)


            # leave only GO_id, filename and pvalues rows in resulting dataframe
            # and convert it to wide format
            resdf <-
              resdf[, c(
                "GO_id",
                "name",
                "namespace",
                "filename",
                "pvalues"
              )]

            resdf <- tidyr::spread(resdf, "filename", "pvalues")

            # move GO_ids from first column to rownames
            rownames(resdf) <- resdf[, 1]
            resdf[, 1] <- NULL
            go_ids <- rownames(resdf)
            go_names <- resdf[, 1]
            go_namespaces <- resdf[, 2]
            resdf <- resdf[, c(-1, -2)]

            # find minimal p-values for each GO term and corresponding interval name
            minp <-
              apply(resdf, 1, function(x)
                c(min(x, na.rm = TRUE), colnames(resdf)[which.min(x)]))

            # from table with genes delete all intervals that are not in minp

            # apply correction on multiple testing
            # and separate GO terms into fold-specific and not fold-specific groups
            # using fdr step 2 threshold
            padj <- stats::p.adjust(minp[1, ], method = p_adjust_method)

            wi_pvaldf <- wi_pvaldf[match(go_ids, wi_pvaldf$GO_id), ]

            result_table <- data.frame(
              ids = go_ids,
              namespace = go_namespaces,
              name = go_names,
              wholeint_pval = wi_pvaldf$pvalue,
              wholeint_padj = wi_pvaldf$padj,
              min_pval = minp[1, ],
              padj = padj,
              interval = minp[2, ],
              stringsAsFactors = FALSE
            )

            object@result_table <- result_table

            return(object)
          }
)
