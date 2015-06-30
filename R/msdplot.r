#' Plot means and standard deviations with a dotplot
#'
#' Function to plot means and confidence limits.
#' @param trait The trait to plot.
#' @param groups The grouping factor.
#' @param data The name of the data frame containing the data.
#' @param conf Probability for the confidence limits or number of standard deviations.
#' @param nmax Maximum number of points for the compulsory dotplot, default is 10.
#' @param dotplot Logical. If \code{TRUE}, a dotplot is shown. If \code{FALSE}
#' it will only suppress the dots if the number of data points is larger than \code{nmax}.
#' @param sort.means Sort for means. Options are \code{"none"}, \code{"increasing"},
#' and \code{"decreasing"}, \code{"none"} by default.
#' @param main Main title.
#' @param xlab Title for x axis.
#' @param ylab Title for y axis.
#' @param colors Color for mean symbols, confidence interval lines, and data points.
#' @param pch Ploting character for means.
#' @param lwd Width for ploting characters for means.
#' @param x.las x axes labels orientation.
#' @param jf Jitter factor for dots.
#' @param dist Horizontal distance between the means and the dots.
#' @author Raul Eyzaguirre
#' @details An alternative to the controversial dynamite plots.
#' If \code{conf} is set to a value greater than or equal to 1, then it is interpreted
#' as number of standard deviations.
#' @return It returns a plot with the means represented by horizontal lines, a vertical
#' line representing a confidence limit or a number of standard deviations,
#' and alternatively the individual data points.
#' @examples
#' # Simulate some data
#' mydata <- data.frame(y = rnorm(50, sample(40:60, 5), sample(5:10, 5)),
#'                      g = rep(1:5, 10))
#'
#' # Draw the plot
#' msdplot("y", "g", mydata)
#' @export
                         
msdplot <- function(trait, groups, data, conf = 0.95, nmax = 10, dotplot = "TRUE",
                    sort.means = "none", main = NULL, xlab = "groups", ylab = "",
                    colors = c("orange", "orange", "black"), pch = 4, lwd = 2,
                    x.las = 1, jf = 0.1, dist = 0.1) {
  
  # Error messages
  
  if (dotplot != "TRUE" & dotplot != "FALSE")
    stop("dotplot argument must be TRUE or FALSE.")    
    
  # Groups as factor
  
    data[, groups] <- as.factor(data[, groups])

  # means and standard deviations

  means <- tapply(data[,trait], data[,groups], mean, na.rm=T)
  sdev <- tapply(data[,trait], data[,groups], sd, na.rm=T)

  resu <- data.frame(means, sdev)

  # compute confidence intervals

  if (conf < 1) {
    resu$n <- tapply(is.na(data[, trait])==0, data[,groups], sum)
    resu$li <- resu$means - qt((1 + conf)/2, resu$n-1) * resu$sdev/sqrt(resu$n)
    resu$ls <- resu$means + qt((1 + conf)/2, resu$n-1) * resu$sdev/sqrt(resu$n)
    msg <- paste("Dotplot with means and ", conf*100, "% confidence limits", sep="")
  } else {
    resu$li <- resu$means - conf * resu$sdev
    resu$ls <- resu$means + conf * resu$sdev
    msg <- paste("Dotplot with means +/-", conf, "standard deviations")
  }
  
  resu$orden <- rownames(resu)

  # sort

  if (sort.means == "increasing")
    resu <- resu[sort(resu$means, index.return=T)$ix,]
  if (sort.means == "decreasing")
    resu <- resu[sort(resu$means, decreasing=T, index.return=T)$ix,]
  if (sort.means %in% c("none", "increasing", "decreasing") == F)
    stop("Invalid value for sort.means")

  # title for plot

  if (is.null(main) == 1)
    main = msg

  # limits for plot
  
  a <- min(resu$li)
  b <- max(resu$ls)
  
  for (i in 1:length(resu$means)){
    subdata <- subset(data, data[, groups] == resu$orden[i])
    if (dotplot == "TRUE" | length(subdata[, trait]) <= nmax){
      a <- min(a, subdata[, trait])
      b <- max(b, subdata[, trait])
    }    
  }

  # draw the plot
  
  plot(seq(1, length(resu$means)), resu$means, xaxt = "n",
       xlab = xlab, ylab = ylab, main = main,
       xlim = c(0.5, length(resu$means) + 0.5), ylim = c(a, b),
       col = colors[1], pch = pch, lwd = lwd)

  axis(1, at = seq(1, length(resu$means)), labels = rownames(resu), las = x.las)

  for (i in 1:length(resu$means)){
    lines(c(i,i), c(resu$li[i], resu$ls[i]), col = colors[2])
    subdata <- subset(data, data[, groups] == resu$orden[i])
    if (dotplot == "TRUE" | length(subdata[, trait]) <= nmax)
      points(jitter(rep(i + dist, length(subdata[, trait])), factor = jf),
             subdata[, trait], col = colors[3])
    }
}
