#' Fold-specificity chart plotting
#'
#' @param FoldSpecTest - object of S4 FoldSpecTest class
#'
#' @export
#'
#' @examples
#' # calculate fold-specificity test for up-regulated genes
#' up_fs <- FoldSpecTest(up_groups)
#' # calculate fold-specificity test for down-regulated genes
#' down_fs <- FoldSpecTest(down_groups)
#' plot(up_fs, down_fs)
setMethod(f = "plot",
          signature = "FoldSpecTest",
          definition = function(x, y, x_text_size = 10){
            up_obj <- x
            down_obj <- y
            wholeintname <- getWholeIntName(up_obj)
            borders <- unlist(strsplit(wholeintname, "-", fixed = TRUE))
            x_labs <- seq(borders[1], borders[2])
            fs_data <- fold_spec_chart_data(up_obj@fstable,
                                            up_obj@nfstable,
                                            down_obj@fstable,
                                            down_obj@nfstable,
                                            wholeintname)
            plot(fold_spec_chart(fs_data, x_labs, x_text_size = x_text_size))
          }
)


#' Create rectangle plot coordinates for each Go term
#'
#' @param intervals - vector of intervals in character representation (e.g. "1-6")
#'
#' @return - matrix contains interval start coordinates in the first column and interval end coordinates
#'           in the second row (note: the start coordinate is deducted by 1 in order to make the appropriate layout for rectangles)
#'
create_intervals_matrix <- function(intervals) {
  output_matrix <- matrix(integer(0), ncol = 2)
  # convert character representation of coordinates to numeric
  for (i in 1:length(intervals)) {
    x <- intervals[i]
    if (grepl("^([1-9][0-9]*)-([1-9][0-9]*)$", x) == TRUE) {
      splitarr <- unlist(strsplit(x, "-", fixed = TRUE))
      start <- as.numeric(splitarr[1]) - 1
      end <- as.numeric(splitarr[2])
      if (start > end)
        stop (paste("start position of interval is greater than end position:", x))
      output_matrix <- rbind(output_matrix, c(start, end))
    } else if (grepl("^([1-9][0-9]*)$", x) == TRUE) {
      start <- as.numeric(x) - 1
      end <- as.numeric(x)
      output_matrix <- rbind(output_matrix, c(start, end))
    } else {
      start <- as.numeric(x)
      end <- as.numeric(x)
      output_matrix <- rbind(output_matrix, c(start, end))
    }
  }
  return(output_matrix)
}

#' Data parser for fold-specificity rectangle plot
#'
#' Create input data for rectangle plot (fold_spec_chart function) using recognize_fs_terms function output as input
#'
#' @param fs_res_up dataframe contains fold-specificty recognition data (recognize_fs_terms function output) for up regulation
#' @param fs_res_down dataframe contains fold-specificty recognition data (recognize_fs_terms function output) for down regulation
#'
#' @return input dataframe for create.fold.spec.chart function
#'
fold_spec_chart_data <-
  function(fs_res_up, fs_res_down, nfs_res_up, nfs_res_down, wholeintname) {

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("tidyr package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    nfs_res_up$interval <- wholeintname
    nfs_res_down$interval <- wholeintname

    nfs_res_up <- nfs_res_up[nfs_res_up$ids %in% fs_res_down$ids, ]
    nfs_res_down <- nfs_res_down[nfs_res_down$ids %in% fs_res_up$ids, ]

    fs_res_up <- rbind(fs_res_up, nfs_res_up)
    fs_res_down <- rbind(fs_res_down, nfs_res_down)

    # add column with regulation type
    fs_res_up$reg <- rep("up", length(fs_res_up[, 1]))
    fs_res_down$reg <- rep("down", length(fs_res_down[, 1]))

    # combine up and down dataframes and
    # spread the resulting one by regulation type
    fs_res_combined <- rbind(fs_res_up, fs_res_down)

    fs_res_combined$name <- paste(fs_res_combined$ids, fs_res_combined$name)

    fs_res_combined <- fs_res_combined[, c(-1, -2, -4, -5, -6, -7)]

    fs_res_combined_spreaded <-
      tidyr::spread(fs_res_combined, "reg", "interval", fill = "0")


    # add zero's if there is no data for certain regulation type
    if (is.null(fs_res_combined_spreaded$up)) {
      fs_res_combined_spreaded$up <-
        rep("0", length(fs_res_combined_spreaded$name))
    }

    if (is.null(fs_res_combined_spreaded$down)) {
      fs_res_combined_spreaded$down <-
        rep("0", length(fs_res_combined_spreaded$name))
    }

    # create matrices with intervals coordinates
    in_mat_up <- FoldGO::create_intervals_matrix(fs_res_combined_spreaded$up)
    in_mat_down <-
      FoldGO::create_intervals_matrix(fs_res_combined_spreaded$down)

    # add start - end coordinates to resulting dataframe
    fs_res_combined_spreaded$start_up <- in_mat_up[, 1]
    fs_res_combined_spreaded$end_up <- in_mat_up[, 2]
    fs_res_combined_spreaded$start_down <- in_mat_down[, 1]
    fs_res_combined_spreaded$end_down <- in_mat_down[, 2]
    fs_res_combined_spreaded <- fs_res_combined_spreaded[, c(-2, -3)]
    fs_res_combined_spreaded <-
      fs_res_combined_spreaded[order(fs_res_combined_spreaded[, 3]), ]

    return(fs_res_combined_spreaded)
  }



#' fold specificity rectangle plot
#'
#' Create rectangle plot for fold specificity data representation
#'
#' @param data dataframe contains GO term names and coordinates for rectangles both for up and down regulation
#' @param interval_labels vector of user defined names for non-overlaping intervals (e.g. c("weak response",...,"strong response"))
#' @param x_text_size size of text for x axis labels
#'
#' @return fold specificity rectangle plot as ggplot object
#' @importFrom ggplot2 ggplot geom_rect scale_x_continuous scale_y_continuous theme geom_hline geom_text coord_flip aes element_blank element_text
#'
fold_spec_chart <- function(data, interval_labels, x_text_size) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # extract y-axis labels
  labels <- as.vector(data[, 1])
  # extract reactangle plot coordinates
  df <- data.frame(c(data[, 2]), c(data[, 3]))
  # create vector of x - axis coordinates
  labs <- c(1:length(data[, 1]))
  # create palindromic vector of y labels
  ylabs <- c(rev(interval_labels), interval_labels)
  limit <- length(interval_labels)
  breaks_limit <- limit - 0.5
  p_up <- ggplot(df, aes(x = labs)) +
    geom_rect(
      aes(
        x = labs,
        xmin = as.numeric(labs) - 0.45,
        xmax = as.numeric(labs) + 0.45,
        ymin = c(data[, 2]),
        ymax = c(data[, 3])
      ),
      fill = "yellow2",
      alpha = 0.8
    ) +
    geom_rect(
      aes(
        x = labs,
        xmin = as.numeric(labs) - 0.45,
        xmax = as.numeric(labs) + 0.45,
        ymin = c(data[, 4]) * -1,
        ymax = c(data[, 5]) * -1
      ),
      fill = "skyblue2",
      alpha = 0.8
    ) +
    scale_x_continuous(breaks = 1:length(labels), labels = labels) +
    scale_y_continuous(
      breaks = seq(-breaks_limit, breaks_limit, 1),
      limits = c(-limit, limit),
      labels = ylabs,
      minor_breaks = seq(-limit, limit, 0.5)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(
        size = x_text_size,
        angle = 45,
        hjust = 1
      )
    ) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               size = 0.2) +
    geom_text(y = 0, aes(label = labels)) +
    coord_flip()
}
