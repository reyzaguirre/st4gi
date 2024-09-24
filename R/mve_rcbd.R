#' Estimation of missing values for a RCBD
#'
#' Function to estimate missing values for a Randomized Complete Block Design (RCBD) by
#' the least squares method.
#' @param y The name of the column for the variable to estimate missing values.
#' @param geno The name of the column that identifies the genotypes.
#' @param rep The name of the column that identifies the replications.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{y}
#' and \code{y.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @examples
#' temp <- met8x12[met8x12$env == "TM80N", ]
#' mve.rcbd("y", "geno", "rep", temp)
#' @export

mve.rcbd <- function(y, geno, rep, dfr, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- ck.rcbd(y, geno, rep, dfr)

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

  y.est <- paste0(y, ".est")
  dfr[, y.est] <- dfr[, y]
  dfr[, "ytemp"] <- dfr[, y]
  mG <- tapply(dfr[, y], dfr[, geno], mean, na.rm = TRUE)
  for (i in 1:length(dfr[, y]))
    if (is.na(dfr[i, y]))
      dfr[i, "ytemp"] <- mG[dfr[i, geno]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(dfr[, y], na.rm = TRUE)
  cont <- 0
  while (cc > max(dfr[, y], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(dfr[, y]))
      if (is.na(dfr[i, y])) {
        dfr[i, "ytemp"] <- dfr[i, y]
        sum1 <- tapply(dfr[, "ytemp"], dfr[, geno], sum, na.rm = TRUE)
        sum2 <- tapply(dfr[, "ytemp"], dfr[, rep], sum, na.rm = TRUE)
        sum3 <- sum(dfr[, "ytemp"], na.rm = TRUE)
        dfr[i, y.est] <- (lc$ng * sum1[dfr[i, geno]] + lc$nrep * sum2[dfr[i, rep]] - sum3) /
          (lc$ng * lc$nrep - lc$ng - lc$nrep + 1)
        dfr[i, "ytemp"] <- dfr[i, y.est]
      }
    lc1 <- lc2
    lc2 <- dfr[is.na(dfr[, y]), y.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  dfr[, c(geno, rep, y, y.est)]
  
}
