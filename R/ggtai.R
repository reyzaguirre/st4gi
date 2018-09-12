#' Tai ggplot
#' 
#' This function produces a Tai's stability analysis plot (Tai, G. C. C., 1971)
#' using the ggplot2 package.
#' @param x An object of class \code{tai}.
#' @param conf Probability for the Tai limits.
#' @details The limits for alpha and lambda are computed using the mean squares from
#' an ANOVA table for a RCBD with blocks nested into environments.See \code{?tai} for
#' additional details.
#' @return It returns the Tai graph for stability analysis using the ggplot2 package.
#' @author Raul Eyzaguirre.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' model.tai <- tai("y", "geno", "env", "rep", met8x12)
#' ggtai(model.tai)
#' @importFrom stats qf qt
#' @export

ggtai <- function(x, conf = 0.95) {
  
  # arguments
  
  trait <- x$Trait
  alpha <- x$Tai_values[, 1]
  lambda <- x$Tai_values[, 2] 
  dat <- as.data.frame(x$Tai_values)
  at <- x$ANOVA
  lc <- x$lc
  
  # plot lambda limits
  
  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, lc$nl[2] - 2,
                           lc$nl[2] * (lc$nl[1] - 1) * (lc$nrep - 1))))
  
  # Prediction interval for alpha
  
  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, lc$nl[2] - 2)
  
  div2 <- (lc$nl[2] - 2) * at[2, 3] - (ta^2 + lc$nl[2] - 2) * at[3, 3]
  
  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments.
            Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (lc$nl[1] - 1) * at[5, 3] * at[2, 3]) / 
                        ((at[2, 3] - at[3, 3]) * div2))^0.5
    amax <- max(c(abs(alpha), pi.alpha))
  }
  
  # Tai ggplot
  
  main <- paste0("Tai stability analysis for ", trait)
  
  gg <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_point(ggplot2::aes(x = lambda, y = alpha)) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(expression(lambda)) +
    ggplot2::ylab(expression(alpha))
  
  # alpha limits
  
  if (div2 > 0) {
    dt2 <- as.data.frame(cbind(lx, pi.alpha))
    gg <- gg +
      ggplot2::geom_path(data = dt2, ggplot2::aes(x = dt2$lx, y = dt2$pi.alpha, col = "black")) +
      ggplot2::geom_path(data = dt2, ggplot2::aes(x = dt2$lx, y = -dt2$pi.alpha, col = "black"))
  }
  
  # lambda limits
  
  gg <- gg +
    ggplot2::geom_vline(xintercept = qf((1 - conf) / 2, lc$nl[2] - 2,
                                        lc$nl[2] * lc$nl[1] * (lc$nrep - 1)),
                        col = "gray") +
    ggplot2::geom_vline(xintercept = qf(1 - (1 - conf) / 2, lc$nl[2] - 2,
                                        lc$nl[2] * lc$nl[1] * (lc$nrep - 1)),
                        col = "gray")
  
  # points
  
  gg <- gg +
    ggrepel::geom_text_repel(ggplot2::aes(x = lambda, y = alpha, label = rownames(dat))) +
    ggplot2::theme(legend.position = 'none')
  
  gg
  
}
