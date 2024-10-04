#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Beep notification
#'
#' Used for internal testing, get a beep notification, e.g. after long run.
#'
#' @param n how many beeps to give in succession
beep <- function(n = 5) {
  for(i in seq(n)){
    system("rundll32 user32.dll, MessageBeep -1")
    Sys.sleep(.5)
  }
  message(paste0("beeped at: ", format(Sys.time(), "%F %H-%M-%S")))
}

#' My custom color palettes.
#'
#' @description This function returns a color palette, based on featured
#' qualitative palletes in base R, many are colorblind friendly.
#' @param type character. Can be 'mixed' for mixed colors based on R
#' qualitative palettes, 'colorblind' for Okabe-Ito palette. For gradients,
#' see ?hcl.colors and notes. Default: mixed
#' @param n integer. How many colors to return.
#' @usage generate_color_palette(type = 'mixed', n = NULL)
#' @note Can also look at creating gradients custom with
#' colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "PRGn"))(100)
#' @return Color palette as character vector.
#' @references https://i.stack.imgur.com/H8syw.png. test
#' https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
#' display.brewer.all(n=50, exact.n=FALSE)
generate_color_palette <- function(type = 'mixed', n = NULL) {

  if (!type %in% c('mixed', 'colorblind')) {
    stop("'type not accepted, please see documentation'")
  }

  if (type == 'mixed') {
    palette <- c("#000000", "#DF536B", "#61D04F", "#2297E6", "#28E2E5",
                          "#CD0BBC", "#F5C710", "#9E9E9E", "#F564E3", "#B79F00",
                          "#E69F00", "#009E73", "#0072B2", "#D55E00", "#666666",
                          "#BEAED4", "#FFFF99", "#F0027F", "#A6761D", "#E31A1C",
                          "#16FF32", "#782AB6", "#E4E1E3", "#E5C494", "#CCCCCC",
                          "#B3CDE3", "#CCEBC5", "#DECBE4", "#FDDAEC", "#FFFFCC")
                          names(palette) <- c("black_R4.1", "salmon_R4.2", "lightgreen_R4.3",
                                              "skyblue_R4.4", "teal_R4.5", "darkpink_R4.6",
                                              "yellow_R4.7", "grey_R4.8", "pink_ggplot2.6",
                                              "gold_ggplot2.7", "orange_OI.2", "darkgreen_OI.4",
                                              "darkblue_OI.6", "darkorange_OI.7", "darkgray_Accent.8",
                                              "violet_Accent.2", "lightyellow_Accent.4",
                                              "brightpink_Accent.6", "brown_Dark2.7", "red_Paired.6",
                                              "brightgeen_Alphabet.7", "purple_Alphabet.4",
                                              "offwhite_Polychrome36.2", "sand_Set2.7",
                                              "lightgrey_Pastel2.8", "pastelblue_Pastel1.2",
                                              "pastelgreen_Pastel1.3", "pastelpurple_Pastel1.4",
                                              "pastelpink_Pastel1.8", "pastelyellow_1.6")
  }
  if (type == 'colorblind') {
    palette <- grDevices::palette.colors(palette = "Okabe-Ito")
  }
  if (!is.null(n)) {
    if (n > length(palette)) {
      message(paste0('WARNING: asked for ', n, ' colors but only ', length(palette), ' are available and selected from the palette'))
      n <- length(palette)
    }
    palette <- palette[1:n]
  }

  # TODO check for diverging color palettes as well: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

  return(palette)
}
