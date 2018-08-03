#' Estimation of missing values for a MET in a RCBD
#'
#' Function to estimate missing values for a Multi Environment Trial (MET) with a
#' Randomized Complete Block Design (RCBD) by the least squares method.
#' @param trait The trait to estimate missing values.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @examples
#' mve.met("y", "geno", "env", "rep", met8x12)
#' @export

mve.met <- function(trait, geno, env, rep, dfr, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- ck.f(trait, c(geno, env), rep, dfr)

  # Error messages

  if (lc$nmis.fac > 0)
    stop("There are missing values for classification factors.")

  if (lc$nt.0 > 0)
    stop("Some GxE cells have zero frequency.")

  if (lc$nrep == 1)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$nt.mult > 0)
    stop("Some genotypes have additional replications.")

  if (lc$nmis == 0)
    stop("There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste0("Too many missing values (", format(lc$pmis * 100, digits = 3), "%)."))

  if (lc$nl[1] < 2 | lc$nl[2] < 2)
    stop("This is not a MET experiment.")

  # Estimation

  trait.est <- paste0(trait, ".est")
  dfr[, trait.est] <- dfr[, trait]
  dfr[, "ytemp"] <- dfr[, trait]
  mGE <- tapply(dfr[, trait], list(dfr[, geno], dfr[, env]), mean, na.rm = TRUE)
  for (i in 1:length(dfr[, trait]))
    if (is.na(dfr[i, trait]))
      dfr[i, "ytemp"] <- mGE[dfr[i, geno], dfr[i, env]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(dfr[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(dfr[, trait], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(dfr[, trait]))
      if (is.na(dfr[i, trait])) {
        dfr[i, "ytemp"] <- dfr[i, trait]
        sum1 <- tapply(dfr[, "ytemp"], list(dfr[, geno], dfr[, env]), sum, na.rm = TRUE)
        sum2 <- tapply(dfr[, "ytemp"], list(dfr[, env], dfr[, rep]), sum, na.rm = TRUE)
        sum3 <- tapply(dfr[, "ytemp"], dfr[, env], sum, na.rm = TRUE)
        dfr[i, trait.est] <- (lc$nl[1] * sum1[dfr[i, geno], dfr[i, env]] +
                                 lc$nrep * sum2[dfr[i, env], dfr[i, rep]] -
                                 sum3[dfr[i, env]]) / (lc$nl[1] * lc$nrep - lc$nl[1] - lc$nrep + 1)
        dfr[i, "ytemp"] <- dfr[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- dfr[is.na(dfr[, trait]), trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  dfr[, c(geno, env, rep, trait, trait.est)]
  
}
