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
  tfreq <- table(subdata[, treat], subdata[, rep])

  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replication
  c3 <- 1 # Check for genotypes with more than one datum in a replication
  c4 <- 1 # Check for missing values. Initial state no missing values
  
  if (min(table(subdata[, treat])) == 0) c1 <- 0 # State 0: there are zeros
  if (nr > 1) c2 <- 1 # State 1: more than one replication
  if (max(tfreq) > 1) c3 <- 0 # State 0: some genotypes with addional data
  if (min(tfreq) == 0) c4 <- 0 # State 0: missing values
  
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
#' the number of replications (\code{nr}), a table with frequencies of valid cases
#' for each combination of the levels of both factors (\code{tfreq}), and a table with
#' frequencies of valid cases for each combination of the levels of both factors in each
#' replication (\code{tfreqr}).
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
  tfreqr <- table(subdata[, A], subdata[, B], subdata[, rep])
  
  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replicate
  c3 <- 1 # Check for genotypes with more than one datum in a replication of one environment
  c4 <- 1 # Check for missing values. Initial state no missing values
    
  if (min(tfreq) == 0) c1 <- 0 # State 0: there are zeros
  if (nr > 1) c2 <- 1 # State 1: more than one replicate
  if (max(tfreqr) > 1) c3 <- 0 # State 0: some genotypes with addional data
  if (min(tfreqr) == 0) c4 <- 0 # State 0: missing values
    
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, nmis = nmis, pmis = pmis,
       na = na, nb = nb, nr = nr, tfreq = tfreq, tfreqr = tfreqr)
}

#' Check rows and columns
#'
#' This function checks that there is only one genotype in each row and column position.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param data The name of the data frame.
#' @return For each replication, a list of row and column positions with more than
#' one genotype.
#' @author Raul Eyzaguirre.
#' @export

check.pos <- function(row, col, rep, data) {
  
  # Number of replications
  
  data[, rep] <- factor(data[, rep])
  lr <- levels(data[, rep])
  nr <- nlevels(data[, rep])
  
  # Check frequencies
  
  for (i in 1:nr) {
    temp <- data[data$rep == lr[i], ]
    ttt <- as.data.frame(table(temp$row, temp$col))
    colnames(ttt) <- c('Row', 'Column', 'Freq')
    dimt <- dim(ttt[ttt$Freq > 1, ])
    cat('------------------------------\n')
    cat('Replication', lr[i], '\n')
    cat('------------------------------\n')
    if (dimt[1] > 0) {
      print(ttt[ttt$Freq > 1, ])
      cat('\n')
    } else {
      cat('OK \n')
      cat('\n')
    }
  }
}

#' Check design
#'
#' Check frequencies for designs with complete replications and one or two factors.
#' @param trait The trait to analyze.
#' @param A Factor A.
#' @param B Factor B.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Information about the balance, missing values, and replications of the design.
#' @author Raul Eyzaguirre.
#' @export

check.design <- function(trait, A, B = NULL, rep, data) {
  
  if (is.null(B)) {
    lc <- check.rcbd(trait, A, rep, data)
    cat('------------------------------\n')
    cat('Check design for trait', trait, '\n')
    cat('------------------------------\n')
    cat('\n')
    if (lc$c1 == 0) {
      cat('There are genotypes without data \n')
      cat('\n')
    }
    if (lc$c2 == 0) {
      cat('There is only one replication \n')
      cat('\n')
    }
    if (lc$c3 == 0) {
      cat('There are genotypes that appear more than once in a given replication \n')
      cat('\n')
    }
    if (lc$c4 == 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '%')
      cat('\n')
    }
    if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1 & lc$c4 == 1) {
      cat('OK \n')
      cat('\n')
    }
  } else {
    lc <- check.2f(trait, A, B, rep, data)
    cat('------------------------------\n')
    cat('Check design for trait', trait, '\n')
    cat('------------------------------\n')
    cat('\n')
    if (lc$c1 == 0) {
      cat('There are genotypes without data in a given environment \n')
      cat('\n')
    }
    if (lc$c2 == 0) {
      cat('There is only one replication \n')
      cat('\n')
    }
    if (lc$c3 == 0) {
      cat('There are genotypes that appear more than once in a given replication \n')
      cat('\n')
    }
    if (lc$c4 == 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '%')
      cat('\n')
    }
    if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1 & lc$c4 == 1) {
      cat('OK \n')
      cat('\n')
    }    
  }
}
