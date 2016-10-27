#' Check data for an ABD
#'
#' This function checks the frequencies of genotypes in an ABD.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return A list of treatments \code{newmat}, a list of checks \code{checks},
#' the number of treatments \code{nt.new}, the number of checks \code{nt.check},
#' the number of missing values for treatments \code{nmis.new}),
#' the number of missing values for checks \code{nmis.check},
#' the number \code{nt.check.0} and the list \code{check.0} of checks without data,
#' the number \code{nt.check.0} and the list \code{check.1} of checks with only one datum,
#' the number of checks with at least two data \code{nt.check.2},
#' and the number of replications \code{nr}.
#' @author Raul Eyzaguirre.
#' @details This function checks the frequencies and number of missing values in an ABD.
#' For an ANOVA in an ABD it is needed at least two checks with at least 2 valid cases each.
#' @export

check.abd <- function(trait, treat, rep, data) {
  
  # Everything as character
  
  data[, treat] <- as.character(data[, treat])
  data[, rep] <- as.character(data[, rep])
  
  # Number of replications
  
  nr <- nlevels(as.factor(data[, rep]))

  # Identify checks and no checks
  
  tfreq <- data.frame(table(data[, treat]))
  checks <- as.character(tfreq[tfreq$Freq > 1, 1])
  newmat <- as.character(tfreq[tfreq$Freq == 1, 1])
  
  # Number of treatments
  
  nt.check <- length(checks)
  nt.new <- length(newmat)
  
  # Number of missing values
  
  temp <- subset(data, !(data[, treat] %in% checks))
  nmis.new <- sum(is.na(temp[, trait]))
  
  temp <- subset(data, data[, treat] %in% checks)
  nmis.check <- sum(is.na(temp[, trait]))
  
  # Number of checks without data, and 1 and more data

  temp[, treat] <- as.factor(temp[, treat])
  temp <- subset(temp, !is.na(temp[, trait]))
  tfreq <- data.frame(table(temp[, treat]))
  
  nt.check.0 <- sum(tfreq$Freq == 0)
  nt.check.1 <- sum(tfreq$Freq == 1)
  nt.check.2 <- sum(tfreq$Freq > 1)
  
  # Checks to remove
  
  check.0 <- NULL
  if (nt.check.0 > 0 )
    check.0 <- as.character(tfreq[tfreq$Freq == 0, 1])
  
  check.1 <- NULL
  if (nt.check.1 > 0)
    check.1 <- as.character(tfreq[tfreq$Freq == 1, 1])

  # Return
  
  list(checks = checks, newmat = newmat, nt.check = nt.check, nt.new = nt.new,
       nmis.check = nmis.check, nmis.new = nmis.new, nr = nr,
       nt.check.0 = nt.check.0, check.0 = check.0,
       nt.check.1 = nt.check.1, check.1 = check.1,
       nt.check.2 = nt.check.2)
}

#' Check data for a RCBD
#'
#' This function checks the frequencies of genotypes in a RCBD.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Four control values (\code{c1}, \code{c2}, \code{c3}, and \code{c4}), the number
#' of missing values \code{nmis}, the proportion of missing values (\code{pmis}), the number
#' of treatments (\code{nt}), the number of replications (\code{nr}), and a table with
#' frequencies of valid cases for each genotype.
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication in a RCBD,
#' if there is any treatment without data or with more data than replications, and
#' if the design is balanced.
#' @export

check.rcbd <- function(trait, treat, rep, data) {
  
  # Everything as factor
  
  data[, treat] <- factor(data[, treat])
  data[, rep] <- factor(data[, rep])
  
  # Number of levels
  
  nt <- nlevels(data[, treat])
  nr <- nlevels(data[, rep])

  # Check frequencies by treat
  
  nmis <- sum(is.na(data[, trait]))
  pmis <- mean(is.na(data[, trait]))
  subdata <- subset(data, is.na(data[, trait]) == 0)
  tfreq <- table(subdata[, treat])

  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replicate
  c3 <- 1 # Check for balance (additional data). Initial state balanced
  c4 <- 1 # Check for missing values. Initial state no missing values
  
  if (min(tfreq) == 0) c1 <- 0 # State 0: there are zeros
  if (max(tfreq) > 1) c2 <- 1 # State 1: more than one replicate
  if (max(tfreq) > nr) c3 <- 0 # State 0: some cells with addional data
  if (min(tfreq) < nr) c4 <- 0 # State 0: missing values
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, nmis = nmis, pmis = pmis,
       nt = nt, nr = nr, tfreq = tfreq)
}

#' Check data for a 2-factor factorial
#'
#' This function checks the frequencies for a 2-factor factorial.
#' @param trait The trait to analyze.
#' @param A Factor A.
#' @param B Factor B.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Four control values (\code{c1}, \code{c2}, \code{c3}, and \code{c4}), the number
#' of missing values \code{nmis}, the proportion of missing values (\code{pmis}), the number
#' of levels of factor A (\code{na}), the number of levels of factor B (\code{nb}),
#' the number of replications (\code{nr}), and a table with frequencies of valid cases
#' for each combination of the levels of both factors.
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication, if there is
#' any combination of the levels of both factors without data or with more data then
#' replications, and if the design is balanced.
#' @export

check.2f <- function(trait, A, B, rep, data) {
  
  # Everything as factor
  
  data[, A] <- factor(data[, A])
  data[, B] <- factor(data[, B])
  data[, rep] <- factor(data[, rep])
  
  # Number of levels
  
  na <- nlevels(data[, A])
  nb <- nlevels(data[, B])
  nr <- nlevels(data[, rep])

  # Check frequencies by A and B
  
  nmis <- sum(is.na(data[, trait]))
  pmis <- mean(is.na(data[, trait]))
  subdata <- subset(data, is.na(data[, trait]) == 0)
  tfreq <- table(subdata[, A], subdata[, B])

  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replicate
  c3 <- 1 # Check for balance (additional data). Initial state balanced
  c4 <- 1 # Check for missing values. Initial state no missing values
    
  if (min(tfreq) == 0) c1 <- 0 # State 0: there are zeros
  if (max(tfreq) > 1) c2 <- 1 # State 1: more than one replicate
  if (max(tfreq) > nr) c3 <- 0 # State 0: some cells with addional data
  if (min(tfreq) < nr) c4 <- 0 # State 0: missing values
    
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, nmis = nmis, pmis = pmis,
       na = na, nb = nb, nr = nr, tfreq = tfreq)
}
