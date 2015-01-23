#' Check data for a RCBD
#'
#' This function checks the frequencies of genotypes in a RCBD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param data The name of the data frame.
#' @return c1, c2, c3, three control values.
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication in a RCBD,
#' if there is any genotype without data, and if the design is balanced.

checkdata01 <- function(trait, geno, data){
  # Check frequencies by geno
  nmis <- sum(is.na(data[,trait]))
  subdata <- subset(data, is.na(data[,trait]) == 0)
  tfreq <- table(subdata[,geno])
  # Controls
  c1 <- 1 # Check for zeros. Initial state no zeros which is good
  c2 <- 0 # Check for replicates. Initial state only one replicate which is bad
  c3 <- 1 # Check for balance. Initial state balanced which is good
  for (i in 1:dim(tfreq)){
    if (tfreq[i] == 0) c1 <- 0 # State 0: there are zeros
    if (tfreq[i] > 1) c2 <- 1 # State 1: more than one replicate
    if (tfreq[i] != tfreq[1]) c3 <- 0 # State 0: unbalanced
  }
  # Return
  list(c1 = c1, c2 = c2, c3 = c3, nmis = nmis)
}

#' Check data for a MET in a RCBD
#'
#' This function checks the frequencies of genotypes in each environment in a RCBD.
#' @param trait The trait to analyze
#' @param geno The genotypes
#' @param env The environments
#' @param data The name of the data frame
#' @return c1, c2, c3, three control values
#' @author Raul Eyzaguirre
#' @details This function checks if there is more than one replication in a RCBD in
#' several environments, if there is any genotype without data for some specific environments,
#' and if the design is balanced.

checkdata02 <- function(trait, geno, env, data){
  # Check frequencies by geno and env
  nmis <- sum(is.na(data[,trait]))
  subdata <- subset(data, is.na(data[,trait]) == 0)
  tfreq <- table(subdata[,geno], subdata[,env])
  # Controls
  c1 <- 1 # Check for zeros. Initial state no zeros which is good
  c2 <- 0 # Check for replicates. Initial state only one replicate which is bad
  c3 <- 1 # Check for balance. Initial state balanced which is good
  for (i in 1:dim(tfreq)[1])
    for (j in 1:dim(tfreq)[2]){
      if (tfreq[i,j] == 0) c1 <- 0 # State 0: there are zeros
      if (tfreq[i,j] > 1) c2 <- 1 # State 1: more than one replicate
      if (tfreq[i,j] != tfreq[1,1]) c3 <- 0 # State 0: unbalanced
    }
  # Return
  list(c1 = c1, c2 = c2, c3 = c3, nmis = nmis)
}
