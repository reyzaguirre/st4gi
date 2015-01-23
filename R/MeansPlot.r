#' Plot means and standard deviations
#'
#' Function to plot means and confidence limits.
#' @param trait The trait to analyze.
#' @param groups The grouping factor.
#' @param data The name of the data frame containing the data.
#' @param conf Probability for the confidence limits or number of standard deviations.
#' @param sort.means Sort for means, \code{"none"} by default.
#' @param main.title Main title.
#' @param x.title Title for x axis.
#' @param y.title Title for y axis.
#' @param col Line color for circles.
#' @param bg Background color for circles.
#' @param col.lines Line color for confidenced interval lines.
#' @author Raul Eyzaguirre
#' @details An alternative to the controversial dynamite plots.
#' If \code{conf} is set to a value greater than or equal to 1, then it is interpreted
#' as number of standard deviations.
#' @return It returns a plot with the means and confidence limits for each group.
#' @examples
#' # Simulate some data
#' mydata <- data.frame(y = rnorm(100, sample(80:120, 10), sample(10:20, 10)), g = rep(1:10, 10))
#'
#' # Draw the plot
#' msdplot("y", "g", mydata)
#' @export

msdplot <- function(trait, groups, data, conf = 0.95, sort.means = "none",
                      main.title = NULL, x.title = "groups", y.title = "",
                      col = "black", bg = "darkorange", col.lines = "black") {

  # means and standard deviations

  means <- tapply(data[,trait], data[,groups], mean, na.rm=T)
  sdev <- tapply(data[,trait], data[,groups], sd, na.rm=T)

  resu <- data.frame(means, sdev)

  # compute confidence intervals

  if (conf < 1) {
    resu$n <- tapply(is.na(data[, trait])==0, data[,groups], sum)
    resu$li <- resu$means - qt((1 + conf)/2, resu$n-1) * resu$sdev/sqrt(resu$n)
    resu$ls <- resu$means + qt((1 + conf)/2, resu$n-1) * resu$sdev/sqrt(resu$n)
    msg <- paste(conf*100, "% confidence limits", sep="")
  } else {
    resu$li <- resu$means - conf * resu$sdev
    resu$ls <- resu$means + conf * resu$sdev
    msg <- paste("Means +/-", conf, "standard deviations")
  }

  # sort

  if (sort.means == "increasing")
    resu <- resu[sort(resu$means, index.return=T)$ix,]
  if (sort.means == "decreasing")
    resu <- resu[sort(resu$means, decreasing=T, index.return=T)$ix,]
  if (sort.means %in% c("none", "increasing", "decreasing") == F)
    stop("Invalid value for sort.means")

  # make plot

  if (is.null(main.title) == 1)
    main.title = msg

  plot(seq(1, length(resu$means)), resu$means, xaxt = "n",
       xlab = x.title, ylab = y.title, main = main.title,
       xlim = c(0.5, length(resu$means) + 0.5), ylim = c(min(resu$li), max(resu$ls)),
       pch = 21, col = col, bg = bg, cex = 2)

  axis(1, at = seq(1, length(resu$means)), labels = rownames(resu), las = 1)

  for (i in 1:length(resu$means)){
    lines(c(i,i), c(resu$li[i], resu$ls[i]), col = col.lines)
  }

}
