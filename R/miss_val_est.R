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
#' temp <- subset(met8x12, env == "TM80N")
#' mve.rcbd("y", "geno", "rep", temp)
#' @export

mve.rcbd <- function(trait, treat, rep, data, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- ck.rcbd(trait, treat, rep, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Some treatments have zero frequency. Remove treatments to proceed.")

  if (lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$c3 == 0)
    stop("Some treatments have additional replications. Remove those replications to proceed.")

  if (lc$c4 == 1)
    stop("There are no missing values to estimate.")

    if (lc$pmis > maxp)
    stop(paste("Too many missing values (",
               format(lc$pmis * 100, digits = 3), "%).", sep = ""))

  # Estimation

  trait.est <- paste(trait, ".est", sep = "")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  mG <- tapply(data[, trait], data[, treat], mean, na.rm = TRUE)
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait]))
      data[i, "ytemp"] <- mG[data[i, treat]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(data[, trait], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(data[, trait]))
      if (is.na(data[i, trait])) {
        data[i, "ytemp"] <- data[i, trait]
        sum1 <- tapply(data[, "ytemp"], data[, treat], sum, na.rm = TRUE)
        sum2 <- tapply(data[, "ytemp"], data[, rep], sum, na.rm = TRUE)
        sum3 <- sum(data[, "ytemp"], na.rm = TRUE)
        data[i, trait.est] <- (lc$nt * sum1[data[i, treat]] + lc$nr * sum2[data[i, rep]] - sum3) /
          (lc$nt * lc$nr - lc$nt - lc$nr + 1)
        data[i, "ytemp"] <- data[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- subset(data, is.na(data[, trait]))[, trait.est]
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
#' str(met8x12)
#'
#' # Estimate the missing values
#' mve.met("y", "geno", "env", "rep", met8x12)
#' @export

mve.met <- function(trait, geno, env, rep, data, maxp = 0.1, tol = 1e-06) {

  # Check data

  lc <- ck.f(trait, c(geno, env), rep, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Some GxE cells have zero frequency.")

  if (lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")

  if (lc$c3 == 0)
    stop("Some genotypes have additional replications.")

  if (lc$c4 == 1)
    stop("There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste("Too many missing values (",
               format(lc$pmis * 100, digits = 3), "%).", sep = ""))

  if (lc$nl[1] < 2 | lc$nl[2] < 2)
    stop("This is not a MET experiment.")

  # Estimation

  trait.est <- paste(trait, ".est", sep = "")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  mGE <- tapply(data[, trait], list(data[, geno], data[, env]), mean, na.rm = TRUE)
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait]))
      data[i, "ytemp"] <- mGE[data[i, geno], data[i, env]]
  lc1 <- array(0, lc$nmis)
  lc2 <- array(0, lc$nmis)
  cc <- max(data[, trait], na.rm = TRUE)
  cont <- 0
  while (cc > max(data[, trait], na.rm = TRUE) * tol & cont < 100) {
    cont <- cont + 1
    for (i in 1:length(data[, trait]))
      if (is.na(data[i, trait])) {
        data[i, "ytemp"] <- data[i, trait]
        sum1 <- tapply(data[, "ytemp"], list(data[, geno], data[, env]), sum, na.rm = TRUE)
        sum2 <- tapply(data[, "ytemp"], list(data[, env], data[, rep]), sum, na.rm = TRUE)
        sum3 <- tapply(data[, "ytemp"], data[, env], sum, na.rm = TRUE)
        data[i, trait.est] <- (lc$nl[1] * sum1[data[i, geno], data[i, env]] +
                                 lc$nr * sum2[data[i, env], data[i, rep]] -
                                 sum3[data[i, env]]) / (lc$nl[1] * lc$nr - lc$nl[1] - lc$nr + 1)
        data[i, "ytemp"] <- data[i, trait.est]
      }
    lc1 <- lc2
    lc2 <- subset(data, is.na(data[, trait]))[, trait.est]
    cc <- max(abs(lc1 - lc2))
  }

  # Return

  data[, c(geno, env, rep, trait, trait.est)]
}

#' Estimation of missing values for a factorial experiment
#'
#' Function to estimate missing values for factorial experiment with a CRD
#' or a RCBD by the least squares method.
#' @param trait The trait to estimate missing values.
#' @param factors The factors.
#' @param rep The replications or blocks.
#' @param design The statistical design, \code{crd} or \code{rcbd}.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @details A \code{data.frame} with data for a factorial experiment with at least two
#' replications and at least one datum for each factor levels' combination must be loaded.
#' Experimental data with only one replication, any factor levels' combination without data,
#' or more missing values than specified in \code{maxp} will generate an error message.
#' @examples
#' # The data
#' str(asc)
#' 
#' # A copy with some missing values
#' temp <- asc
#' temp$dm[c(3, 10, 115)] <- NA
#'
#' # Estimate the missing values
#' mve.f("dm", c("geno", "treat"), "rep", "crd", temp)
#' @export

mve.f <- function(trait, factors, rep, design = c("crd", "rcbd"), data, maxp = 0.1, tol = 1e-06) {
  
  # match arguments
  
  design <- match.arg(design)
  
  # Check data
  
  lc <- ck.f(trait, factors, rep, data)
  
  # Error messages
  
  if (lc$c1 == 0)
    stop("Some factor levels' combinations have zero frequency.")
  
  if (lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")
  
  if (lc$c3 == 0)
    stop("Some factor levels' combinations have additional replications.")

  if (lc$c4 == 1)
    stop("The data set is balanced. There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste("Too many missing values (",
               format(lc$pmis * 100, digits = 3), "%).", sep = ""))
  
  if (sum(lc$nl < 2) > 0)
    stop("There are factors with only one level. This is is not a factorial experiment.")
  
  # Number of factors
  
  nf <- length(factors)

  # Get a copy of trait for estimated values and for temporary values
  
  trait.est <- paste(trait, ".est", sep = "")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  
  # Create expression for list of factors
  
  lf.expr <- 'list(data[, factors[1]]'
  
  for (i in 2:nf)
    lf.expr <- paste(lf.expr, ', data[, factors[', i, ']]', sep = "")
  
  lf.expr <- paste(lf.expr, ')', sep = "")
  
  # Compute means over replications

  tmeans <- tapply(data[, trait], eval(parse(text = lf.expr)), mean, na.rm = TRUE)

  # Store means in ytemp for missing values
  
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait])) {
      expr <- paste('tmeans[data[', i, ', factors[1]]', sep = "")
      for (j in 2:nf)
        expr <- paste(expr, ', data[', i, ', factors[', j, ']]', sep = "")
      expr <- paste(expr, ']', sep = "")
      data[i, "ytemp"] <- eval(parse(text = expr))
    }
  
  # Estimate missing values for a crd
  
  if (design == "crd")
    data[, trait.est] <- data[, "ytemp"]
  
  # Estimate missing values for a rcbd
  
  if (design == "rcbd"){
    
    # Total number of treatments
    
    nt <- prod(lc$nl)
    
    # Vectors of estimated missing values to check convergence
    
    emv1 <- array(0, lc$nmis)
    emv2 <- array(0, lc$nmis)
    
    # Value to control convergence
    
    cc <- max(data[, trait], na.rm = TRUE)
    
    # Maximum number of iteration control
    
    cont <- 0
    
    # Iterate
    
    while (cc > max(data[, trait], na.rm = TRUE) * tol & cont < 100) {
      
      cont <- cont + 1
      
      for (i in 1:length(data[, trait]))
        if (is.na(data[i, trait])) {
      
          # Compute sums
          
          data[i, "ytemp"] <- NA
          sum1 <- tapply(data[, "ytemp"], eval(parse(text = lf.expr)), sum, na.rm = TRUE)
          sum2 <- tapply(data[, "ytemp"], data[, rep], sum, na.rm = TRUE)
          sum3 <- sum(data[, "ytemp"], na.rm = TRUE)
          
          # Get estimate
          
          expr <- paste('sum1[data[', i, ', factors[1]]', sep = "")
          for (j in 2:nf)
            expr <- paste(expr, ', data[', i, ', factors[', j, ']]', sep = "")
          expr <- paste(expr, ']', sep = "")
          
          mv.num <- nt * eval(parse(text = expr)) + lc$nr * sum2[data[i, rep]] - sum3
          mv.den <- nt * lc$nr - nt - lc$nr + 1

          data[i, trait.est] <- mv.num / mv.den
          
          data[i, "ytemp"] <- data[i, trait.est]
        }
      
      emv1 <- emv2
      emv2 <- subset(data, is.na(data[, trait]) == 1)[, trait.est]
      cc <- max(abs(emv1 - emv2))
    }
  }
  
  # Return
  
  data[, c(factors, rep, trait, trait.est)]
}
