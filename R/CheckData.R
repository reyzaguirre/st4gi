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

checkdata01 <- function(trait, treat, rep, data) {
  
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

checkdata02 <- function(trait, geno, env, rep, data) {
  
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
