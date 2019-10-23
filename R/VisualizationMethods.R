#' Fold-change specific GO Profile chart plotting
#'
#' @param x - object of S4 FoldSpecTest class with up-regulated genes
#' @param y - object of S4 FoldSpecTest class with down-regulated genes
#' @param x_text_size - x axis labels size
#' @return - Fold-change specific GO Profile plot
#' @export
#'
#' @examples
#' # calculate fold-specificity test for up-regulated genes
#' up_fs <- FoldSpecTest(up_annotobj)
#' # calculate fold-specificity test for down-regulated genes
#' down_fs <- FoldSpecTest(down_annotobj)
#' plot(up_fs, down_fs)
setMethod(f = "plot",
          signature = "FoldSpecTest",
          definition = function(x, y, x_text_size = 10, legend_ratio = c(1,5)){
            up_obj <- x
            down_obj <- y
            wholeintname <- getWholeIntName(up_obj)
            #borders <- unlist(strsplit(wholeintname, "-", fixed = TRUE))
            #if (is.null(x_labs)) {
            #  x_labs <- seq(borders[1], borders[2])
            #}
            fs_data <- fold_spec_chart_data(up_obj@fstable,
                                            down_obj@fstable,
                                            up_obj@nfstable,
                                            down_obj@nfstable,
                                            wholeintname)
            fs_chart <- fold_spec_chart(fs_data, up_obj@fold_borders,
                                        down_obj@fold_borders, x_text_size = x_text_size)
            leg_chart <- legend_plot(as.numeric(sub("1-", "", wholeintname)))
            p <- gridExtra::grid.arrange(leg_chart, fs_chart, heights = legend_ratio)
            plot(p)
          }
)


# Create rectangle plot coordinates for each Go term
##############################PARAMS#########################################
# intervals - vector of intervals in character representation (e.g. "1-6")
#
# return - matrix contains interval start coordinates in the first column and
#          interval end coordinates in the second row (note: the start coordinate
#         is deducted by 1 in order to make the appropriate layout for rectangles)
#
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

#  Data parser for fold-specificity rectangle plot
#
#  Create input data for rectangle plot (fold_spec_chart function) using recognize_fs_terms function output as input
######################################PARAMS#################################################################
#  param fs_res_up - dataframe contains fold-specific GO terms and related data for up-regulation
#  param fs_res_down - dataframe contains fold-specific GO terms and related data for down-regulation
#  param nfs_res_up - dataframe contains not fold-specific GO terms and related data for up-regulation
#  param nfs_res_down - dataframe contains not fold-specific GO terms and related data for down-regulation
#  param wholeintname - name of the interval containing all differentially expressed genes
#  input dataframe for create.fold.spec.chart function
fold_spec_chart_data <-
  function(fs_res_up, fs_res_down, nfs_res_up, nfs_res_down, wholeintname) {

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("tidyr package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    #nfs_res_up$interval <- wholeintname
    #nfs_res_down$interval <- wholeintname

    #nfs_res_up <- nfs_res_up[nfs_res_up$ids %in% fs_res_down$ids, ]
    #nfs_res_down <- nfs_res_down[nfs_res_down$ids %in% fs_res_up$ids, ]

    #fs_res_up <- rbind(fs_res_up, nfs_res_up)
    #fs_res_down <- rbind(fs_res_down, nfs_res_down)

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
    in_mat_up <-
      create_intervals_matrix(fs_res_combined_spreaded$up)
    in_mat_down <-
      create_intervals_matrix(fs_res_combined_spreaded$down)

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



# fold specificity rectangle plot
#
# Create rectangle plot for fold specificity data representation
###################################PARAMS####################################################################################
#  data dataframe contains GO term names and coordinates for rectangles both for up and down regulation
#  param interval_labels vector of user defined names for non-overlaping intervals (e.g. c("weak response",...,"strong response"))
#  param x_text_size size of text for x axis labels
#
# return fold specificity rectangle plot as ggplot object
fold_spec_chart <- function(data, up_borders, down_borders, x_text_size = 10) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  up_borders <- formatC(up_borders[-1], digits = 3)
  down_borders <- formatC(down_borders[-1], digits = 3)
  interval_labels <- c(rev(down_borders), "0.00", up_borders)
  # extract y-axis labels
  labels <- as.vector(data[, 1])
  # extract reactangle plot coordinates
  df <- data.frame(c(data[, 2]), c(data[, 3]))
  # create vector of x - axis coordinates
  labs <- c(1:length(data[, 1]))
  # create palindromic vector of y labels
  #ylabs <- interval_labels
  limit <- (length(interval_labels) - 1)/2
  breaks_limit <- limit
  p_up <- ggplot(df, aes(x = labs)) +
    geom_hline(yintercept = c(c(-limit:-1), c(1:limit)),
               linetype = "solid",
               colour = "white",
               size = 0.2) +
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
      labels = interval_labels,
      minor_breaks = NULL
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
    coord_flip() + labs(y = "Log fold change")
}

# Fold quantiles legend plot
#
# param n - number of quantiles
#
# return - ggplot2 object
legend_plot <- function(n) {
  if (n < 11) {
    num_labs <- seq(0.5, n - 0.5, 1)
    p <- ggplot() + geom_polygon(
      data = data.frame(x = c(0, n, n),
                        y = c(0, 0, 1)),
      mapping = aes(x = x, y = y),
      fill = "yellow2",
      alpha = 0.8
    ) +
      geom_polygon(
        data = data.frame(x = c(0, -n, -n),
                          y = c(0, 0, 1)),
        mapping = aes(x = x, y = y),
        fill = "skyblue2",
        alpha = 0.8
      ) +
      geom_vline(xintercept = c(-n:n),
                 linetype = "dashed",
                 size = 0.2) +
      geom_text(x = c(-n, 0, n), y = -0.01,
                aes(label = c("strongest", "weakest", "strongest"))) +
      geom_text(
        x = c(-num_labs, num_labs),
        y = 0.55,
        aes(label = as.character(c(1:n, 1:n))),
        size = 5
      ) +
      geom_label(x = c(-n / 2, n / 2), y = 0.25,
                 aes(label = c("down regulation", "up regulation"))) +
      geom_text(x = 0,
                y = 0.75,
                aes(label = "Quantiles"),
                size = 7) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()
      )
    return(p)
  } else {
    num_labs <- c(seq(0.5, 2.5, 1), seq(7.5, 9.5, 1))
    p <- ggplot() + geom_polygon(
      data = data.frame(x = c(0, 10, 10),
                        y = c(0, 0, 1)),
      mapping = aes(x = x, y = y),
      fill = "yellow2",
      alpha = 0.8
    ) +
      geom_polygon(
        data = data.frame(x = c(0, -10, -10),
                          y = c(0, 0, 1)),
        mapping = aes(x = x, y = y),
        fill = "skyblue2",
        alpha = 0.8
      ) +
      geom_vline(xintercept = c(-10:-7, -3:3, 7:10),
                 linetype = "dashed",
                 size = 0.2) +
      geom_text(x = c(-10, 0, 10), y = -0.01,
                aes(label = c("strongest", "weakest", "strongest"))) +
      geom_text(
        x = c(-num_labs, num_labs),
        y = 0.55,
        aes(label = as.character(c(1:3, (n-2):n, 1:3, (n-2):n))),
        size = 5
      ) +
      geom_label(x = c(-5, 5), y = 0.25,
                 aes(label = c("down regulation", "up regulation"))) +
      geom_text(x = 0,
                y = 0.75,
                aes(label = "Quantiles"),
                size = 7) +
      geom_text(x = -5,
                y = 0.55,
                aes(label = "..."),
                size = 7) +
      geom_text(x = 5,
                y = 0.55,
                aes(label = "..."),
                size = 7) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()
      )
    return(p)
  }
}
