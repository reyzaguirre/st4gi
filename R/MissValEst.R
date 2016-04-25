#' Estimation of missing values for a RCBD
#'
#' Function to estimate missing values for a Randomized Complete Block Design (RCBD) by
#' the least squares method.
#' @param trait The trait to estimate missing values.
#' @param treat The treatments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @details A \code{data.frame} with data for a RCBD with at least two replications
#' and at least one datum for each treatment must be loaded. Experimental data
#' with only one replication, any treatment without data, or more missing values than
#' specified in \code{maxp} will generate an error message.
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Choose one environment
#' temp <- subset(met8x12, env == "TM80N")
#'
#' # Missing value in the first row
#' head(temp)
#'
#' # Estimate the missing value
#' mveb("y", "geno", "rep", temp)
#' @export

mveb <- function(trait, treat, rep, data, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- check.rcbd(trait, treat, rep, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Some treatments have zero frequency. Remove treatments to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1)
    stop("The data set is balanced. There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste("Too many missing values (",
               format(lc$pmis * 100, digits = 3), "%).", sep = ""))

  # Estimation

  trait.est <- paste(trait, ".est", sep = "")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  mG <- tapply(data[, trait], data[, treat], mean, na.rm = TRUE)
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait]) == 1) data[i, "ytemp"] <- mG[data[i, treat]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(data[, trait], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(data[, trait]))
      if (is.na(data[i, trait]) == 1) {
        data[i, "ytemp"] <- data[i, trait]
        sum1 <- tapply(data[, "ytemp"], data[, treat], sum, na.rm = TRUE)
        sum2 <- tapply(data[, "ytemp"], data[, rep], sum, na.rm = TRUE)
        sum3 <- sum(data[, "ytemp"], na.rm = TRUE)
        data[i, trait.est] <- (lc$nt * sum1[data[i, treat]] + lc$nr * sum2[data[i, rep]] - sum3) /
          (lc$nt * lc$nr - lc$nt - lc$nr + 1)
        data[i, "ytemp"] <- data[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- subset(data, is.na(data[, trait]) == 1)[, trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  data[, c(treat, rep, trait, trait.est)]
}

#' Estimation of missing values for a MET in a RCBD
#'
#' Function to estimate missing values for a Multi Environment Trial (MET) with a
#' Randomized Complete Block Design (RCBD) by the least squares method.
#' @param trait The trait to estimate missing values.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @details A \code{data.frame} with data for a MET in a RCBD with at least two replications
#' and at least one datum for each treatment must be loaded. Experimental data
#' with only one replication, any treatment without data, or more missing values than
#' specified in \code{maxp} will generate an error message.
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Estimate the missing values
#' mvemet("y", "geno", "env", "rep", met8x12)
#' @export

mvemet <- function(trait, geno, env, rep, data, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- check.met(trait, geno, env, rep, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1)
    stop("The data set is balanced. There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste("Too many missing values (",
               format(lc$pmis * 100, digits = 3), "%).", sep = ""))

  if (lc$ng < 2 | lc$ne < 2)
    stop("This is not a MET experiment.")

  # Estimation

  trait.est <- paste(trait, ".est", sep = "")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  mGE <- tapply(data[, trait], list(data[, geno], data[, env]), mean, na.rm = TRUE)
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait]) == 1) data[i, "ytemp"] <- mGE[data[i, geno], data[i, env]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(data[, trait], na.rm = TRUE) * tol & cont<100) {
    cont <- cont + 1
    for (i in 1:length(data[, trait]))
      if (is.na(data[i, trait]) == 1) {
        data[i, "ytemp"] <- data[i, trait]
        sum1 <- tapply(data[, "ytemp"], list(data[, geno], data[, env]), sum, na.rm = TRUE)
        sum2 <- tapply(data[, "ytemp"], list(data[, env], data[, rep]), sum, na.rm = TRUE)
        sum3 <- tapply(data[, "ytemp"], data[, env], sum, na.rm = TRUE)
        data[i, trait.est] <- (lc$ng * sum1[data[i, geno], data[i, env]] +
                                 lc$nr * sum2[data[i, env], data[i, rep]] -
                                 sum3[data[i, env]]) / (lc$ng * lc$nr - lc$ng - lc$nr + 1)
        data[i, "ytemp"] <- data[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- subset(data, is.na(data[, trait]) == 1)[, trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  data[, c(geno, env, rep, trait, trait.est)]
}
