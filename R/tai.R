#' Tai's stability analysis
#'
#' This function runs Tai's stability analysis (Tai, G. C. C., 1971).
#' It assumes a RCBD with fixed effects for genotypes and random effects for environments.
#' @param dfr The name of the data frame.
#' @param y The name of the column for the variable to analyze.
#' @param geno The name of the column that identifies the genotypes.
#' @param env The name of the column that identifies the environments.
#' @param rep The name of the column that identifies the replications.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @details If the data set is unbalanced, a warning is produced.
#' @return It returns the alpha and lambda values for each genotype for the Tai
#' stability analysis.
#' @author Raul Eyzaguirre.
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
#' @examples
#' model.tai <- tai(met8x12, "y", "geno", "env", "rep")
#' model.tai$Tai_values
#' @export

tai <- function(dfr, y, geno, env, rep, maxp = 0.1) {

  # Everything as character

  dfr[, geno] <- as.character(dfr[, geno])
  dfr[, env] <- as.character(dfr[, env])
  dfr[, rep] <- as.character(dfr[, rep])

  # Check data

  lc <- ck.f(dfr, y, c(geno, env), rep)

  # Error messages and warnings

  if (lc$nl[1] < 3 | lc$nl[2] < 3)
    stop("You need at least 3 genotypes and 3 environments to run Tai")

  # Estimate missing values and report errors from mve.met
  
  y.est <- paste0(y, ".est")
  
  if (lc$nt.0 > 0 | lc$nrep == 1 | lc$nt.mult > 0 | lc$nmis > 0 | lc$nmis.fac > 0) {
    dfr[, y] <- mve.met(dfr, y, geno, env, rep, maxp, tol = 1e-06)[, y.est]
    warning(paste0("The data set is unbalanced, ",
                   format(lc$pmis * 100, digits = 3),
                   "% missing values estimated."))
  }

  # Compute interaction effects matrix

  int.mean <- tapply(dfr[, y], list(dfr[, geno], dfr[, env]), mean, na.rm = TRUE)

  overall.mean <- mean(int.mean)
  env.mean <- apply(int.mean, 2, mean)
  geno.mean <- apply(int.mean, 1, mean)
  int.eff <- int.mean + overall.mean
  int.eff <- int.eff - geno.mean
  int.eff <- t(t(int.eff)- env.mean)

  # ANOVA

  model <- aov(dfr[, y] ~ dfr[, geno] + dfr[, env] +
                 dfr[, rep] %in% dfr[, env] + dfr[, geno]:dfr[, env])
  at <- anova(model)
  
  # Correction for missing values if any
  
  if (lc$nmis > 0) {
    at[5, 1] <- at[5, 1] - lc$nmis
    at[5, 3] <- at[5, 2] / at[5, 1]
  }
  
  # Compute Tai values alpha and lambda

  slgl <- int.eff
  slgl <- t(t(slgl) * (env.mean - overall.mean) / (lc$nl[2] - 1))
  alpha <- apply(slgl, 1, sum) / (at[2, 3] - at[3, 3]) * lc$nl[1] * lc$nrep

  s2gl <- int.eff
  s2gl <- s2gl^2 / (lc$nl[2] - 1)
  lambda <- (apply(s2gl, 1, sum) - alpha * apply(slgl, 1, sum)) / 
    (lc$nl[1] - 1) / at[5, 3] * lc$nl[1] * lc$nrep
  lambda[lambda < 0] <- 0

  # Output

  output <- list(Variable = y, Tai_values = cbind(alpha, lambda), ANOVA = at, lc = lc)
  
  class(output) <- "st4gi_tai"
  invisible(output)
  
}
