#' Tai's stability analysis
#'
#' This function runs Tai's stability analysis (Tai, G. C. C., 1971).
#' It assumes a RCBD with fixed effects for genotypes and random effects for environments.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @author Raul Eyzaguirre.
#' @details If the data set is unbalanced, a warning is produced.
#' @return It returns the alpha and lambda values for each genotype for the Tai
#' stability analysis.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' model.tai <- tai("y", "geno", "env", "rep", met8x12)
#' model.tai$Tai_values
#' @export

tai <- function(trait, geno, env, rep, data, maxp = 0.1) {

  # Everything as factor

  data[, geno] <- factor(data[, geno])
  data[, env] <- factor(data[, env])
  data[, rep] <- factor(data[, rep])

  # Check data

  lc <- check.AxB(trait, geno, env, rep, data)

  # Error messages and warnings

  if (lc$c1 == 0)
    stop("Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$na < 2 | lc$nb < 2)
    stop("This is not a MET experiment.")

  if (lc$na < 3 | lc$nb < 3)
    stop("You need at least 3 genotypes and 3 environments to run Tai")

  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 0) {
    data[, trait] <- mvemet(trait, geno, env, rep, data, maxp, tol = 1e-06)[, 5]
    warning(paste("The data set is unbalanced, ",
                  format(lc$pmis * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }

  # Compute interaction effects matrix

  int.mean <- tapply(data[, trait], list(data[, geno], data[, env]), mean, na.rm = TRUE)

  overall.mean <- mean(int.mean)
  env.mean <- apply(int.mean, 2, mean)
  geno.mean <- apply(int.mean, 1, mean)
  int.eff <- int.mean + overall.mean
  int.eff <- int.eff - geno.mean
  int.eff <- t(t(int.eff)- env.mean)

  # ANOVA

  model <- aov(data[, trait] ~ data[, geno] + data[, env] +
                 data[, rep] %in% data[, env] + data[, geno]:data[, env])
  at <- anova(model)
  
  # Correction for missing values if any
  
  if (lc$nmis > 0) {
    at[5, 1] <- at[5, 1] - lc$nmis
    at[5, 3] <- at[5, 2] / at[5, 1]
  }
  
  # Compute Tai values alpha and lambda

  slgl <- int.eff
  slgl <- t(t(slgl) * (env.mean - overall.mean) / (lc$nb - 1))
  alpha <- apply(slgl, 1, sum) / (at[2, 3] - at[3, 3]) * lc$na * lc$nr

  s2gl <- int.eff
  s2gl <- s2gl^2 / (lc$nb - 1)
  lambda <- (apply(s2gl, 1, sum) - alpha * apply(slgl, 1, sum)) / 
    (lc$na - 1) / at[5, 3] * lc$na * lc$nr
  lambda[lambda < 0] <- 0

  # Output

  output <- list(Trait = trait, Tai_values = cbind(alpha, lambda), ANOVA = at, lc = lc)
  
  class(output) <- "tai"
  invisible(output)
}

#' Tai plot
#' 
#' This function produces a Tai's stability analysis plot (Tai, G. C. C., 1971).
#' @param x An object of class \code{tai}.
#' @param conf Probability for the Tai limits.
#' @param title Main title for plot.
#' @param color Color for symbols, labels and lines.
#' @param size Relative size for symbols and labels.
#' @param ... Additional plot arguments.
#' @author Raul Eyzaguirre.
#' @details The limits for alpha and lambda are computed using the mean squares from
#' an ANOVA table for a RCBD with blocks nested into environments.See \code{?tai} for
#' additional details.
#' @return It returns the Tai graph for stability analysis.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' model.tai <- tai("y", "geno", "env", "rep", met8x12)
#' plot(model.tai)
#' @export

plot.tai <- function(x, conf = 0.95, title = NULL, color = c("darkorange", "black", "gray"),
                     size = c(1, 1), ...) {

  # arguments
  
  trait <- x$Trait
  alpha <- x$Tai_values[, 1]
  lambda <- x$Tai_values[, 2] 
  at <- x$ANOVA
  lc <- x$lc

  # plot lambda limits
  
  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, lc$nb - 2,
                           lc$nb * (lc$na - 1) * (lc$nr - 1)))) * 1.1
  
  # Prediction interval for alpha
  
  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, lc$nb - 2)
  
  div2 <- (lc$nb - 2) * at[2, 3] - (ta^2 + lc$nb - 2) * at[3, 3]
  
  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments. Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (lc$na - 1) * at[5, 3] * at[2, 3]) / ((at[2, 3] - at[3, 3]) * div2))^0.5
    amax <- max(c(abs(alpha), pi.alpha))
  }

  # Tai plot
  
  if (is.null(title))
    title <- paste("Tai stability analysis for ", trait, sep = "")
  
  plot(1, type = "n", xlim = c(-0.05 * lmax, lmax), ylim = c(-amax, amax),
       main = title, xlab = expression(lambda), ylab = expression(alpha))
  points(lambda, alpha, col = color[1], lwd = 2, pch = 4, cex = size[1])
  text(lambda, alpha, labels = names(alpha), col = color[2], pos = 1,
       offset = 0.3, cex = size[2])
  if (div2 > 0) {
    points(lx, pi.alpha, type = "l", lty = 5, col = color[3])
    points(lx, -pi.alpha, type = "l", lty = 5, col = color[3])
  }
  abline(v = qf((1 - conf) / 2, lc$nb - 2, lc$nb * lc$na * (lc$nr - 1)),
         lty = 5, col = color[3])
  abline(v = qf(1 - (1 - conf) / 2, lc$nb - 2, lc$nb * lc$na * (lc$nr - 1)),
         lty = 5, col = color[3])
}

#' Tai ggplot
#' 
#' This function produces a Tai's stability analysis plot (Tai, G. C. C., 1971)
#' using the ggplot2 package.
#' @param x An object of class \code{tai}.
#' @param conf Probability for the Tai limits.
#' @param title Main title for plot.
#' @author Raul Eyzaguirre.
#' @details The limits for alpha and lambda are computed using the mean squares from
#' an ANOVA table for a RCBD with blocks nested into environments.See \code{?tai} for
#' additional details.
#' @return It returns the Tai graph for stability analysis using the ggplot2 package.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' model.tai <- tai("y", "geno", "env", "rep", met8x12)
#' ggtai(model.tai)
#' @importFrom stats qf qt
#' @export

ggtai <- function(x, conf = 0.95, title = NULL) {
  
  # arguments
  
  trait <- x$Trait
  alpha <- x$Tai_values[, 1]
  lambda <- x$Tai_values[, 2] 
  dat <- as.data.frame(x$Tai_values)
  dat$geno <- rownames(dat)
  at <- x$ANOVA
  lc <- x$lc
  
  # plot lambda limits
  
  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, lc$nb - 2,
                           lc$nb * (lc$na - 1) * (lc$nr - 1)))) * 1.1
  
  # Prediction interval for alpha
  
  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, lc$nb - 2)
  
  div2 <- (lc$nb - 2) * at[2, 3] - (ta^2 + lc$nb - 2) * at[3, 3]
  
  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments. Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (lc$na - 1) * at[5, 3] * at[2, 3]) / ((at[2, 3] - at[3, 3]) * div2))^0.5
    amax <- max(c(abs(alpha), pi.alpha))
  }
  
  # Tai ggplot
  
  if (is.null(title))
    title <- paste("Tai stability analysis for ", trait, sep = "")
  
  gg <- ggplot2::ggplot(data = dat, ggplot2::aes(x = lambda, y = alpha)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(expression(lambda)) +
    ggplot2::ylab(expression(alpha)) +
    ggplot2::coord_cartesian(xlim = c(-0.05 * lmax, lmax), ylim = c(-amax, amax))
  
  # alpha limits
  
  if (div2 > 0) {
    dt2 <- as.data.frame(cbind(lx, pi.alpha))
    gg <- gg +
      ggplot2::geom_path(data = dt2, ggplot2::aes(x = dt2$lx, y = dt2$pi.alpha, col = "black")) +
      ggplot2::geom_path(data = dt2, ggplot2::aes(x = dt2$lx, y = -dt2$pi.alpha, col = "black"))
  }
  
  # lambda limits
  
  gg <- gg +
    ggplot2::geom_vline(xintercept = qf((1 - conf) / 2, lc$nb - 2, lc$nb * lc$na * (lc$nr - 1)),
                        col = "gray") +
    ggplot2::geom_vline(xintercept = qf(1 - (1 - conf) / 2, lc$nb - 2, lc$nb * lc$na * (lc$nr - 1)),
                        col = "gray")
  
  # points
  
  gg <- gg +
    ggrepel::geom_text_repel(ggplot2::aes(label = dat$geno, col = "blue")) +
    ggplot2::theme(legend.position = 'none')
  
  gg
}
