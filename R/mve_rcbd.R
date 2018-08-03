#' Estimation of missing values for a RCBD
#'
#' Function to estimate missing values for a Randomized Complete Block Design (RCBD) by
#' the least squares method.
#' @param trait The trait to estimate missing values.
#' @param geno The genotypes.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @examples
#' temp <- met8x12[met8x12$env == "TM80N", ]
#' mve.rcbd("y", "geno", "rep", temp)
#' @export

mve.rcbd <- function(trait, geno, rep, dfr, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- ck.rcbd(trait, geno, rep, dfr)

  # Error messages

  if (lc$nmis.fac > 0)
    stop("There are missing values for classification factors.")

  if (lc$ng.0 > 0)
    stop("Some genotypes have zero frequency.")

  if (lc$nrep == 1)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$ng.mult > 0)
    stop("Some genotypes have additional replications.")

  if (lc$nmis == 0)
    stop("There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste0("Too many missing values (", format(lc$pmis * 100, digits = 3), "%)."))

  # Estimation

  trait.est <- paste0(trait, ".est")
  dfr[, trait.est] <- dfr[, trait]
  dfr[, "ytemp"] <- dfr[, trait]
  mG <- tapply(dfr[, trait], dfr[, geno], mean, na.rm = TRUE)
  for (i in 1:length(dfr[, trait]))
    if (is.na(dfr[i, trait]))
      dfr[i, "ytemp"] <- mG[dfr[i, geno]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(dfr[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(dfr[, trait], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(dfr[, trait]))
      if (is.na(dfr[i, trait])) {
        dfr[i, "ytemp"] <- dfr[i, trait]
        sum1 <- tapply(dfr[, "ytemp"], dfr[, geno], sum, na.rm = TRUE)
        sum2 <- tapply(dfr[, "ytemp"], dfr[, rep], sum, na.rm = TRUE)
        sum3 <- sum(dfr[, "ytemp"], na.rm = TRUE)
        dfr[i, trait.est] <- (lc$ng * sum1[dfr[i, geno]] + lc$nrep * sum2[dfr[i, rep]] - sum3) /
          (lc$ng * lc$nrep - lc$ng - lc$nrep + 1)
        dfr[i, "ytemp"] <- dfr[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- dfr[is.na(dfr[, trait]), trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  dfr[, c(geno, rep, trait, trait.est)]
  
}
