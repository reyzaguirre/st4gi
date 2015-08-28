#' Tai's stability analysis
#'
#' This function runs Tai's stability analysis (Tai, G. C. C., 1971).
#' It assumes a RCBD with fixed effects for genotypes and random effects for environments.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param conf Probability for the Tai limits.
#' @param title Main title for plot.
#' @param color Color for symbols, labels and lines.
#' @param ... Additional graphic parameters.
#' @author Raul Eyzaguirre
#' @details The limits for alpha and lambda are computed using the mean squares from
#' an ANOVA table for a RCBD with blocks nested into environments. If the data set is
#' unbalanced, a warning is produced.
#' @return It returns the Tai graph for stability analysis and the values of alpha
#' and lambda for each genotype.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Run Tai for trait y
#' tai("y", "geno", "env", "rep", met8x12)
#' @export

tai <- function(trait, geno, env, rep, data, maxp = 0.1, conf = 0.95,
                title = NULL, color = c("darkorange", "black", "gray"), ...) {

  # Everything as factor

  data[, geno] <- factor(data[, geno])
  data[, env] <- factor(data[, env])
  data[, rep] <- factor(data[, rep])

  # Check data

  lc <- checkdata02(trait, geno, env, data)

  # Error messages and warnings

  geno.num <- nlevels(data[, geno])
  env.num <- nlevels(data[, env])
  rep.num <- lc$rep.num

  if (lc$c1 == 0)
    stop("Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (geno.num < 2 | env.num < 2)
    stop("This is not a MET experiment.")

  if (geno.num < 3 | env.num < 3)
    stop("You need at least 3 genotypes and 3 environments to run Tai")

  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 0) {
    est.data <- mvemet(trait, geno, env, rep, data, maxp, tol = 1e-06)
    data[, trait] <- est.data$new.data[, 5]
    nmis <- est.data$est.num
    warning(paste("The data set is unbalanced, ",
                  format(est.data$est.prop * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  } else {
    nmis <- 0
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
  
  if (nmis > 0) {
    at[5, 1] <- at[5, 1] - nmis
    at[5, 3] <- at[5, 2] / at[5, 1]
  }
  
  # Compute Tai values alpha and lambda

  slgl <- int.eff
  slgl <- t(t(slgl) * (env.mean - overall.mean) / (env.num - 1))
  alpha <- apply(slgl, 1, sum) / (at[2, 3] - at[3, 3]) * geno.num * rep.num

  s2gl <- int.eff
  s2gl <- s2gl^2 / (env.num - 1)
  lambda <- (apply(s2gl, 1, sum) - alpha * apply(slgl, 1, sum)) / 
    (geno.num - 1) / at[5, 3] * geno.num * rep.num
  lambda[lambda < 0] <- 0

  # plot lambda limits

  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, env.num - 2, env.num * (geno.num - 1) * (rep.num - 1)))) * 1.1

  # Prediction interval for alpha

  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, env.num - 2)
  
  div2 <- (env.num - 2) * at[2, 3] - (ta^2 + env.num - 2) * at[3, 3]

  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments. Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (geno.num - 1) * at[5, 3] * at[2, 3]) / ((at[2, 3] - at[3, 3]) * div2))^.5
    amax <- max(c(abs(alpha), pi.alpha))
  }
  
  # Tai plot

  if (is.null(title))
    title <- paste("Tai stability analysis for ", trait, sep = "")

  plot(1, type = "n", xlim = c(-0.05 * lmax, lmax), ylim = c(-amax, amax),
       main = title, xlab = expression(lambda), ylab = expression(alpha), ...)
  points(lambda, alpha, col = color[1], lwd = 2, pch = 4, ...)
  text(lambda, alpha, labels = names(alpha), col = color[2], pos = 1, offset = .3)
  if (div2 > 0) {
    points(lx, pi.alpha, type = "l", lty = 3, col = color[3])
    points(lx, -pi.alpha, type = "l", lty = 3, col = color[3])
  }
  abline(v = qf((1 - conf) / 2, env.num - 2, env.num * geno.num * (rep.num - 1)),
         lty = 3, col = color[3])
  abline(v = qf(1 - (1 - conf) / 2, env.num - 2, env.num * geno.num * (rep.num - 1)),
         lty = 3, col = color[3])

  # Output

  coords <- cbind(alpha, lambda)
  return(coords)

}
