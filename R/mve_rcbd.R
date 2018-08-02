#' Estimation of missing values for a RCBD
#'
#' Function to estimate missing values for a Randomized Complete Block Design (RCBD) by
#' the least squares method.
#' @param trait The trait to estimate missing values.
#' @param geno The genotypes.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @details A \code{data.frame} with data for a RCBD with at least two replications
#' and at least one datum for each genotype must be loaded. Experimental data
#' with only one replication, any genotype without data, or more missing values than
#' specified in \code{maxp} will generate an error message.
#' @examples
#' temp <- met8x12[met8x12$env == "TM80N", ]
#' mve.rcbd("y", "geno", "rep", temp)
#' @export

mve.rcbd <- function(trait, geno, rep, data, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- ck.rcbd(trait, geno, rep, data)

  # Error messages

  if (lc$nmis.fact > 0)
    stop("There are missing values for classification factors.")

  if (lc$c1 == 0)
    stop("Some genotypes have zero frequency.")

  if (lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$c3 == 0)
    stop("Some genotypes have additional replications.")

  if (lc$c4 == 1)
    stop("There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste0("Too many missing values (", format(lc$pmis * 100, digits = 3), "%)."))

  # Estimation

  trait.est <- paste0(trait, ".est")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  mG <- tapply(data[, trait], data[, geno], mean, na.rm = TRUE)
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait]))
      data[i, "ytemp"] <- mG[data[i, geno]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(data[, trait], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(data[, trait]))
      if (is.na(data[i, trait])) {
        data[i, "ytemp"] <- data[i, trait]
        sum1 <- tapply(data[, "ytemp"], data[, geno], sum, na.rm = TRUE)
        sum2 <- tapply(data[, "ytemp"], data[, rep], sum, na.rm = TRUE)
        sum3 <- sum(data[, "ytemp"], na.rm = TRUE)
        data[i, trait.est] <- (lc$ng * sum1[data[i, geno]] + lc$nr * sum2[data[i, rep]] - sum3) /
          (lc$ng * lc$nr - lc$ng - lc$nr + 1)
        data[i, "ytemp"] <- data[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- data[is.na(data[, trait]), trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  data[, c(geno, rep, trait, trait.est)]
  
}
