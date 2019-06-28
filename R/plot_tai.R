#' Tai plot
#' 
#' This function produces a Tai's stability analysis plot (Tai, G. C. C., 1971)
#' using the base system and ggplot2.
#' @param x An object of class \code{tai}.
#' @param conf Probability for the Tai limits.
#' @param graph.type \code{"base"} or \code{"ggplot"}.
#' @param color Color for symbols, labels and lines (Only for the base system plot).
#' @param ... Additional plot arguments (Only for the base system plot).
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
#' plot(model.tai, graph.type = "ggplot")
#' @importFrom stats qf qt
#' @export

plot.st4gi_tai <- function(x, conf = 0.95, graph.type = c("base", "ggplot"),
                           color = c("darkorange", "black", "gray"), ...) {

  # match arguments
  
  graph.type <- match.arg(graph.type)

  # arguments
  
  trait <- x$Trait
  alpha <- x$Tai_values[, 1]
  lambda <- x$Tai_values[, 2] 
  at <- x$ANOVA
  lc <- x$lc
  
  if (graph.type == "ggplot")
    dat <- as.data.frame(x$Tai_values)
  
  # plot lambda limits
  
  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, lc$nl[2] - 2,
                           lc$nl[2] * (lc$nl[1] - 1) * (lc$nrep - 1))))
  
  if (graph.type == "base")
    lmax <- lmax * 1.1
  
  # Prediction interval for alpha
  
  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, lc$nl[2] - 2)
  
  div2 <- (lc$nl[2] - 2) * at[2, 3] - (ta^2 + lc$nl[2] - 2) * at[3, 3]
  
  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments.
            Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (lc$nl[1] - 1) * at[5, 3] * at[2, 3]) / ((at[2, 3] - at[3, 3]) * div2))^0.5
    amax <- max(c(abs(alpha), pi.alpha))
  }

  # Tai plot
  
  main <- paste0("Tai stability analysis for ", trait)
  
  if (graph.type == "base") {
    
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
    
  } else {
    
    gg <- ggplot2::ggplot(data = dat) +
      ggplot2::geom_point(ggplot2::aes(x = lambda, y = alpha)) +
      ggplot2::ggtitle(main) +
      ggplot2::xlab(expression(lambda)) +
      ggplot2::ylab(expression(alpha))
    
    # Alpha limits
    
    if (div2 > 0) {
      dt2 <- as.data.frame(cbind(lx, pi.alpha))
      gg <- gg +
        ggplot2::geom_path(data = dt2, ggplot2::aes(x = dt2$lx, y = dt2$pi.alpha, col = "black")) +
        ggplot2::geom_path(data = dt2, ggplot2::aes(x = dt2$lx, y = -dt2$pi.alpha, col = "black"))
    }
    
    # Lambda limits
    
    gg <- gg +
      ggplot2::geom_vline(xintercept = qf((1 - conf) / 2, lc$nl[2] - 2,
                                          lc$nl[2] * lc$nl[1] * (lc$nrep - 1)),
                          col = "gray") +
      ggplot2::geom_vline(xintercept = qf(1 - (1 - conf) / 2, lc$nl[2] - 2,
                                          lc$nl[2] * lc$nl[1] * (lc$nrep - 1)),
                          col = "gray")
    
    # Points
    
    gg <- gg +
      ggrepel::geom_text_repel(ggplot2::aes(x = lambda, y = alpha, label = rownames(dat))) +
      ggplot2::theme(legend.position = 'none')
    
    gg
    
  }
  
}
