#' Check data for an ABD
#'
#' This function checks the frequencies of genotypes in an ABD.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return A list of treatments \code{newmat}, a list of checks \code{checks},
#' the number of treatments \code{nt.new}, the number of checks \code{nt.check},
#' the number of missing values for treatments \code{nmis.new}), the number of missing
#' values for checks \code{nmis.check}, the number of checks without data \code{nt.check.0},
#' the number of checks with only one datum \code{nt.check.1}, the number of checks with
#' at least two data \code{nt.check.2}, and the number of replications \code{nr}.
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
  checks <- as.character(subset(tfreq, Freq > 1)[, 1])
  newmat <- as.character(subset(tfreq, Freq == 1)[, 1])
  
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
    
  # Return
  
  list(checks = checks, newmat = newmat, nt.check = nt.check, nt.new = nt.new,
       nmis.check = nmis.check, nmis.new = nmis.new, nr = nr,
       nt.check.0 = nt.check.0, nt.check.1 = nt.check.1, nt.check.2 = nt.check.2)
}

#' Check data for a RCBD
#'
#' This function checks the frequencies of genotypes in a RCBD.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Three control values (\code{c1}, \code{c2}, and \code{c3}), the number of
#' missing values \code{nmis}, the proportion of missing values (\code{pmis}), the number
#' of treatments (\code{nt}), and the number of replications (\code{nr}).
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication in a RCBD,
#' if there is any treatment without data, and if the design is balanced.
#' @export

check.rcbd <- function(trait, treat, rep, data) {
  
  # Everything as factor
  
  data[, treat] <- factor(data[, treat])
  data[, rep] <- factor(data[, rep])
  
  nt <- nlevels(data[, treat])
  nr <- nlevels(data[, rep])
  
  # Check frequencies by treat
  
  nmis <- sum(is.na(data[, trait]))
  pmis <- mean(is.na(data[, trait]))
  subdata <- subset(data, is.na(data[, trait]) == 0)
  tfreq <- table(subdata[, treat])

  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros which is good
  c2 <- 0 # Check for replicates. Initial state only one replicate which is bad
  c3 <- 1 # Check for balance. Initial state balanced which is good
  
  if (min(tfreq) == 0) c1 <- 0 # State 0: there are zeros
  if (max(tfreq) > 1) c2 <- 1 # State 1: more than one replicate
  if (min(tfreq) != max(tfreq)) c3 <- 0 # State 0: unbalanced
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, nmis = nmis, pmis = pmis, nt = nt, nr = nr)
}

#' Check data for a MET in a RCBD
#'
#' This function checks the frequencies of genotypes in each environment in a RCBD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Three control values (\code{c1}, \code{c2}, and \code{c3}), the number of
#' missing values \code{nmis}, the proportion of missing values (\code{pmis}), the number
#' of genotypes (\code{ng}), the number of environments (\code{ne}), and the number of
#' replications (\code{nr}).
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication in a RCBD in
#' several environments, if there is any genotype without data for some specific environments,
#' and if the design is balanced.
#' @export

check.met <- function(trait, geno, env, rep, data) {
  
  # Everything as factor
  
  data[, geno] <- factor(data[, geno])
  data[, env] <- factor(data[, env])
  data[, rep] <- factor(data[, rep])
  
  ng <- nlevels(data[, geno])
  ne <- nlevels(data[, env])
  nr <- nlevels(data[, rep])

  # Check frequencies by geno and env
  
  nmis <- sum(is.na(data[, trait]))
  pmis <- mean(is.na(data[, trait]))
  subdata <- subset(data, is.na(data[, trait]) == 0)
  tfreq <- table(subdata[, geno], subdata[, env])

  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros which is good
  c2 <- 0 # Check for replicates. Initial state only one replicate which is bad
  c3 <- 1 # Check for balance. Initial state balanced which is good
  
  if (min(tfreq) == 0) c1 <- 0 # State 0: there are zeros
  if (max(tfreq) > 1) c2 <- 1 # State 1: more than one replicate
  if (min(tfreq) != max(tfreq)) c3 <- 0 # State 0: unbalanced
    
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, nmis = nmis, pmis = pmis, ng = ng, ne = ne, nr = nr)
}
