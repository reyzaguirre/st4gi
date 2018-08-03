#' Tai plot
#' 
#' This function produces a Tai's stability analysis plot (Tai, G. C. C., 1971).
#' @param x An object of class \code{tai}.
#' @param conf Probability for the Tai limits.
#' @param color Color for symbols, labels and lines.
#' @param ... Additional plot arguments.
#' @details The limits for alpha and lambda are computed using the mean squares from
#' an ANOVA table for a RCBD with blocks nested into environments.See \code{?tai} for
#' additional details.
#' @return It returns the Tai graph for stability analysis.
#' @author Raul Eyzaguirre.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' model.tai <- tai("y", "geno", "env", "rep", met8x12)
#' plot(model.tai)
#' @export

plot.tai <- function(x, conf = 0.95, color = c("darkorange", "black", "gray"), ...) {

  # arguments
  
  trait <- x$Trait
  alpha <- x$Tai_values[, 1]
  lambda <- x$Tai_values[, 2] 
  at <- x$ANOVA
  lc <- x$lc

  # plot lambda limits
  
  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, lc$nl[2] - 2,
                           lc$nl[2] * (lc$nl[1] - 1) * (lc$nrep - 1)))) * 1.1
  
  # Prediction interval for alpha
  
  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, lc$nl[2] - 2)
  
  div2 <- (lc$nl[2] - 2) * at[2, 3] - (ta^2 + lc$nl[2] - 2) * at[3, 3]
  
  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments. Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (lc$nl[1] - 1) * at[5, 3] * at[2, 3]) / ((at[2, 3] - at[3, 3]) * div2))^0.5
    amax <- max(c(abs(alpha), pi.alpha))
  }

  # Tai plot
  
  main <- paste0("Tai stability analysis for ", trait)
  
  plot(1, type = "n", xlim = c(-0.05 * lmax, lmax), ylim = c(-amax, amax),
       main = main, xlab = expression(lambda), ylab = expression(alpha))
  points(lambda, alpha, col = color[1], lwd = 2, pch = 4)
  text(lambda, alpha, labels = names(alpha), col = color[2], pos = 1, offset = 0.3)
  if (div2 > 0) {
    points(lx, pi.alpha, type = "l", lty = 5, col = color[3])
    points(lx, -pi.alpha, type = "l", lty = 5, col = color[3])
  }
  abline(v = qf((1 - conf) / 2, lc$nl[2] - 2, lc$nl[2] * lc$nl[1] * (lc$nrep - 1)),
         lty = 5, col = color[3])
  abline(v = qf(1 - (1 - conf) / 2, lc$nl[2] - 2, lc$nl[2] * lc$nl[1] * (lc$nrep - 1)),
         lty = 5, col = color[3])
  
}
