#' Check rows and columns
#'
#' This function checks that there is only one genotype in each row and column position.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param data The name of the data frame.
#' @return For each replication, the number (\code{nplot}) and list (\code{lplot})
#' of plots (unique row and column position) with more than one genotype.
#' @author Raul Eyzaguirre.
#' @export

check.rc <- function(row, col, rep = NULL, data) {
  
  # Number of replications
  
  if (is.null(rep)) {
    data[, "rep"] <- 1
    rep <- "rep"
  }
  data[, rep] <- factor(data[, rep])
  lr <- levels(data[, rep])
  nr <- nlevels(data[, rep])
  
  # Check row and column
  
  nplot <- NULL
  lplot <- list()
  
  for (i in 1:nr) {
    
    # Compute frequencies
    
    temp <- data[data[, rep] == lr[i], ]
    ttt <- as.data.frame(table(temp[, row], temp[, col]))
    colnames(ttt) <- c('Row', 'Column', 'Freq')
    
    # Save number of plots with problems
    
    nplot[i] <- dim(ttt[ttt$Freq > 1, ])[1]
    
    # Save list of plots with problems if any
    
    if (nplot[i] > 0)
      lplot[[i]] <- ttt[ttt$Freq > 1, ] 
  }
  
  # Return
  
  list(nplot = nplot, lplot = lplot, lr = lr, nr = nr)
}

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
  temp[, treat] <- as.factor(temp[, treat])
  temp[, rep] <- as.factor(temp[, rep])
  temp <- subset(temp, !is.na(temp[, trait]))
  tfreq <- table(temp[, treat], temp[, rep])
  nmis.check <- sum(tfreq == 0)

  # Number of checks without data, with 1 and more data

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
  
  subdata <- subset(data, !is.na(data[, trait]))
  tfreq <- table(subdata[, treat], subdata[, rep])
  nmis <- sum(tfreq == 0)
  pmis <- mean(tfreq == 0)
  
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
  subdata <- subset(data, !is.na(data[, trait]))
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

#' Check data for a Wescott design
#'
#' This function checks the grid of checks on the Wescott design and
#' the number of missing values.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param data The name of the data frame.
#' @return Five control values (\code{c1}, \code{c2}, \code{c3}, \code{c4}, \code{c5})
#' for the grid of checks, the number of missing values for checks (\code{nmis.check})
#' and genotypes \code{nmis}, and the proportion of missing values for checks
#' (\code{pmis.check}) and genotypes (\code{pmis}).
#' @author Raul Eyzaguirre.
#' @details This function checks the grid of checks for the Wescoot design anc
#' calculates the number of missing values.
#' @export

check.wd <- function(trait, geno, ch1, ch2, row, col, ncb, data) {
  
  # Checks
  
  checks <- c(ch1, ch2)
  
  # Numbers and characters
  
  data[, row] <- as.numeric(data[, row])
  data[, col] <- as.numeric(data[, col])
  data[, geno] <- as.character(data[, geno])
  
  # Number of rows and columns
  
  nr.min <- min(data[, row])
  nc.min <- min(data[, col])
  
  nc.max <- max(data[, col])
  
  # Controls
  
  c1 <- 0 # All column checks with checks
  c2 <- 0 # Last column with checks
  c3 <- 0 # All column genotypes with genotypes
  c4 <- 0 # Alternating checks without in correlative row order
  c5 <- 0 # All genotypes with checks to the left and right

  # Columns with checks
  
  cch <- seq(nc.min, nc.max, ncb + 1)
  
  # Check columns with checks
  
  if (sum(!(data[data[, col] %in% cch, geno] %in% checks)) > 0)
    c1 <- 1
  
  # Last column with checks
  
  if (max(cch) != nc.max)
    c2 <- 1
  
  # Check columns with genotypes
  
  if (sum(data[!(data[, col] %in% cch), geno] %in% checks) > 0)
    c3 <- 1
  
  # Alternating checks
  
  for (i in cch)
    for (j in (min(data[data[, col] == i, row]) + 1):max(data[data[, col] == i, row]))
      if (data[data[, col] == i & data[, row] == j, geno] == data[data[, col] == i & data[, row] == j - 1, geno])
        c4 <- 1
  
  for (i in 2:length(cch))
    if (data[data[, col] == cch[i] & data[, row] == nr.min, geno] == data[data[, col] == cch[i - 1] & data[, row] == nr.min, geno])
      c4 <- 1
  
  # All genotypes must have one check to the left and one to the right
  
  for(i in 1:dim(data)[1]) {
    rows <- data[i, row]
    columns <- (data[i, col] - ncb):(data[i, col] + ncb)
    temp <- data[data[, row] == rows & data[, col] %in% columns, geno]
    if (sum(temp %in% checks) == 0)
      c5 <- 1
  }

  # Number of missing values for checks
  
  temp <- data[data[, col] %in% cch, ]
  total.check <- dim(temp)[1]
  temp <- temp[is.na(temp[, trait]), ]
  nmis.check <- dim(temp)[1]
  pmis.check <- nmis.check / total.check
  
  # Number of missing values for genotypes
  
  temp <- data[!(data[, col] %in% cch), ]
  total.geno <- dim(temp)[1]
  temp <- temp[is.na(temp[, trait]), ]
  nmis <- dim(temp)[1]
  pmis <- nmis / total.geno
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, nmis = nmis, pmis = pmis,
       nmis.check = nmis.check, pmis.check = pmis.check)
}
