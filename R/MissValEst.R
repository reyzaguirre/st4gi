#' Estimation of missing values for a RCBD
#'
#' Function to estimate missing values for a Randomized Complete Block Design (RCBD) by
#' the least squares method.
#' @param trait The trait to estimate missing values.
#' @param geno The genotypes.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 5\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with name \code{new.data} and the number
#' \code{est.num} and proportion \code{est.prop} of estimated missing values.
#' The \code{new.data} data frame contains the experimental layout and columns
#' \code{trait} and \code{trait.est} with the original data and the original data plus the
#' estimated values.
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
#' temp <- subset(met8x12, env=="TM80N")
#'
#' # Missing value in the first row
#' head(temp)
#'
#' # Estimate the missing value
#' mveb("y", "geno", "rep", temp)
#' @export

mveb <- function(trait, geno, rep, data, maxp = 0.05, tol = 1e-06){

  # Everything as factor

  data[,geno] <- factor(data[,geno])
  data[,rep] <- factor(data[,rep])

  # Check data

  lc <- checkdata01(trait, geno, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Error: Some genotypes have zero frequency. Remove genotypes to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    stop("Error: There is only one replication. Inference is not possible with one replication.")

  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1)
    stop("Error: The data set is balanced. There are no missing values to estimate.")

  est.p <- mean(is.na(data[,trait]))
  if (est.p > maxp)
    stop(paste("Error: Too many missing values (",
               format(est.p*100, digits = 3), "%).", sep = ""))

  # Estimation

  G <- nlevels(data[,geno])
  R <- nlevels(data[,rep])

  trait.est <- paste(trait, ".est", sep = "")
  data[,trait.est] <- data[,trait]
  data[,"ytemp"] <- data[,trait]
  mG <- tapply(data[,trait], data[,geno], mean, na.rm = T)
  for (i in 1:length(data[,trait]))
    if (is.na(data[i,trait]) == 1) data[i,"ytemp"] <- mG[data[i,geno]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[,trait], na.rm = T)
  cont <- 0
  while (cc > max(data[,trait], na.rm = T)*tol & cont<100){
    cont <- cont+1
    for (i in 1:length(data[,trait]))
      if (is.na(data[i,trait]) == 1){
        data[i,"ytemp"] <- data[i,trait]
        sum1 <- tapply(data[,"ytemp"], data[,geno], sum, na.rm = T)
        sum2 <- tapply(data[,"ytemp"], data[,rep], sum, na.rm = T)
        sum3 <- sum(data[,"ytemp"], na.rm = T)
        data[i,trait.est] <- (G*sum1[data[i,geno]] + R*sum2[data[i,rep]] - sum3)/
                             (G*R - G - R + 1)
        data[i,"ytemp"] <- data[i,trait.est]
      }
    lc1 <- lc2
    lc2 <- subset(data, is.na(data[,trait]) == 1)[,trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  list(new.data = data[,c(geno,rep,trait,trait.est)],
       est.num = lc$nmis, est.prop = est.p)
}

#' Estimation of missing values for a MET in a RCBD
#'
#' Function to estimate missing values for a Multi Environment Trial (MET) with a
#' Randomized Complete Block Design (RCBD) by the least squares method.
#' @param trait The trait to estimate missing values.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 5\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with name \code{new.data} and the number
#' \code{est.num} and proportion \code{est.prop} of estimated missing values.
#' The \code{new.data} data frame contains the experimental layout and columns
#' \code{trait} and \code{trait.est} with the original data and the original data plus the
#' estimated values.
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

mvemet <- function(trait, geno, env, rep, data, maxp = 0.05, tol = 1e-06){

  # Everything as factor

  data[,geno] <- factor(data[,geno])
  data[,env] <- factor(data[,env])
  data[,rep] <- factor(data[,rep])

  # Check data

  lc <- checkdata02(trait, geno, env, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Error: Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    stop("Error: There is only one replication. Inference is not possible with one replication.")

  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1)
    stop("Error: The data set is balanced. There are no missing values to estimate.")

  est.p <- mean(is.na(data[,trait]))
  if (est.p > maxp)
    stop(paste("Error: Too many missing values (",
               format(est.p*100, digits = 3), "%).", sep = ""))

  G <- nlevels(data[,geno])
  E <- nlevels(data[,env])
  R <- nlevels(data[,rep])

  if (G < 2 | E < 2)
    stop(paste("Error: This is not a MET experiment."))

  # Estimation

  trait.est <- paste(trait, ".est", sep = "")
  data[,trait.est] <- data[,trait]
  data[,"ytemp"] <- data[,trait]
  mGE <- tapply(data[,trait], list(data[,geno], data[,env]), mean, na.rm = T)
  for (i in 1:length(data[,trait]))
    if (is.na(data[i,trait]) == 1) data[i,"ytemp"] <- mGE[data[i,geno], data[i,env]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[,trait], na.rm = T)
  cont <- 0
  while (cc > max(data[,trait], na.rm = T)*tol & cont<100){
    cont <- cont+1
    for (i in 1:length(data[,trait]))
      if (is.na(data[i,trait]) == 1){
        data[i,"ytemp"] <- data[i,trait]
        sum1 <- tapply(data[,"ytemp"], list(data[,geno], data[,env]), sum, na.rm = T)
        sum2 <- tapply(data[,"ytemp"], list(data[,env], data[,rep]), sum, na.rm = T)
        sum3 <- tapply(data[,"ytemp"], data[,env], sum, na.rm = T)
        data[i,trait.est] <- (G*sum1[data[i,geno], data[i,env]] +
                                R*sum2[data[i,env], data[i,rep]] -
                                sum3[data[i,env]]) / (G*R - G - R + 1)
        data[i,"ytemp"] <- data[i,trait.est]
      }
    lc1 <- lc2
    lc2 <- subset(data, is.na(data[,trait]) == 1)[,trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  list(new.data = data[,c(geno,env,rep,trait,trait.est)],
       est.num = lc$nmis, est.prop = est.p)
}
