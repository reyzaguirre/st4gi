#' Plot means and standard deviations with a dotplot
#'
#' Function to plot means and confidence limits.
#' @param trait The trait to plot.
#' @param groups The grouping factor.
#' @param dfr The name of the data frame.
#' @param conf Probability for the confidence limits or number of standard deviations.
#' @param dotplot Logical. If \code{TRUE}, a dotplot is shown.
#' @param sort.means Sort for means. Options are \code{"none"}, \code{"increasing"},
#' and \code{"decreasing"}, \code{"none"} by default.
#' @param main.title Main title.
#' @param color Color for mean symbols, confidence interval lines, and data points.
#' @param x.las x axes labels orientation.
#' @param jf Jitter factor for dots.
#' @param hsep Horizontal separation between the means and the dots.
#' @param ... Additional plot arguments.
#' @details An alternative to the controversial dynamite plots.
#' If \code{conf} is set to a value greater than or equal to 1, then it is interpreted
#' as number of standard deviations.
#' @return It returns a plot with the means represented by horizontal lines, a vertical
#' line representing a confidence limit or a number of standard deviations,
#' and alternatively the individual data points.
#' @author Raul Eyzaguirre
#' @examples
#' # Simulate some data
#' dfr <- data.frame(y = rnorm(50, sample(40:60, 5), sample(5:10, 5)),
#'                   g = rep(1:5, 10))
#' msdplot("y", "g", dfr, lwd = 2, pch = 4)
#' @importFrom graphics axis lines plot points
#' @export
                         
msdplot <- function(trait, groups, dfr, conf = 0.95, dotplot = TRUE,
                    sort.means = c("none", "increasing", "decreasing"),
                    main.title = NULL, color = c("orange", "orange", "black"),
                    x.las = 1, jf = 0.1, hsep = 0.1, ...) {

  # Match arguments

  sort.means <- match.arg(sort.means)
  
  # Delete missing values

  dfr <- dfr[!is.na(dfr[, trait]), ]
  
  # Means and standard deviations

  means <- tapply(dfr[, trait], dfr[, groups], mean, na.rm = TRUE)
  sdev <- tapply(dfr[, trait], dfr[, groups], sd, na.rm = TRUE)

  resu <- data.frame(means, sdev)

  # Compute confidence intervals

  if (conf < 1) {
    resu$n <- table(dfr[, groups])
    resu$li <- resu$means - qt((1 + conf) / 2, resu$n - 1) * resu$sdev / sqrt(resu$n)
    resu$ls <- resu$means + qt((1 + conf) / 2, resu$n - 1) * resu$sdev / sqrt(resu$n)
    msg <- paste0("Dotplot with means and ", conf * 100, "% confidence limits")
  } else {
    resu$li <- resu$means - conf * resu$sdev
    resu$ls <- resu$means + conf * resu$sdev
    msg <- paste("Dotplot with means +/-", conf, "standard deviations")
  }
  
  if (!is.null(main.title))
    msg <- main.title
  
  resu$orden <- rownames(resu)

  # sort

  if (sort.means == "increasing")
    resu <- resu[sort(resu$means, index.return = TRUE)$ix, ]
  if (sort.means == "decreasing")
    resu <- resu[sort(resu$means, decreasing = TRUE, index.return = TRUE)$ix, ]
  if (sort.means %in% c("none", "increasing", "decreasing") == FALSE)
    stop("Invalid value for sort.means")

  # limits for plot
  
  a <- min(resu$li, na.rm = T)
  b <- max(resu$ls, na.rm = T)
  
  for (i in 1:length(resu$means)) {
    subdata <- dfr[dfr[, groups] == resu$orden[i], ]
    if (dotplot == "TRUE") {
      a <- min(a, subdata[, trait])
      b <- max(b, subdata[, trait])
    }    
  }

  # draw the plot
  
  group <- seq(1, length(resu$means))
  y <- resu$means
  
  plot(group, y, xaxt = "n", main = msg, xlim = c(0.5, length(resu$means) + 0.5),
       ylim = c(a, b), col = color[1], ...)

  axis(1, at = seq(1, length(resu$means)), labels = rownames(resu), las = x.las)

  for (i in 1:length(resu$means)) {
    lines(c(i, i), c(resu$li[i], resu$ls[i]), col = color[2])
    subdata <- dfr[dfr[, groups] == resu$orden[i], ]
    if (dotplot == "TRUE")
      points(jitter(rep(i + hsep, length(subdata[, trait])), factor = jf),
             subdata[, trait], col = color[3])
  }
  
}
