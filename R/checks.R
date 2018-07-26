#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param data The name of the data frame.
#' @return For each replication, the number (\code{nplot}) and list (\code{lplot})
#' of plots (unique row and column position) with more than one genotype, the 
#' number (\code{nr}) and list (\code{lr}) of replications, and the number of rows
#' in the data frame with missing values for factors (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @export

ck.pos <- function(row, col, rep = NULL, data) {
  
  # Define replications
  
  if (is.null(rep)) {
    data[, "rep"] <- 1
    rep <- "rep"
  }
  
  # Check and delete rows with missing values for factors
  
  nmis.fact <- nrow(data[is.na(data[, row]) | is.na(data[, col]) | is.na(data[, rep]), ])
  
  if (nmis.fact > 0)
    data <- data[!(is.na(data[, row]) | is.na(data[, col]) | is.na(data[, rep])), ]

  # Number of replications and levels

  lr <- sort(unique(data[, rep]))
  nr <- length(lr)
  
  # Number and list of plots with more than one genotype
  
  nplot <- NULL
  lplot <- list()
  
  # Check row and column

  for (i in 1:nr) {
    
    # Compute frequencies
    
    temp <- data[data[, rep] == lr[i], ]
    ttt <- as.data.frame(table(temp[, row], temp[, col]))
    colnames(ttt) <- c('Row', 'Column', 'Freq')
    
    # Number of plots with problems
    
    nplot[i] <- dim(ttt[ttt$Freq > 1, ])[1]
    
    # List of plots with problems if any
    
    if (nplot[i] > 0)
      lplot[[i]] <- ttt[ttt$Freq > 1, ] 
  }
  
  # Return
  
  list(nplot = nplot, lplot = lplot, nr = nr, lr = lr, nmis.fact = nmis.fact)
  
}

#' Check data for an ABD
#'
#' This function checks the frequencies of genotypes in an ABD.
#' @param trait The trait to analyze.
#' @param geno The genotypes including checks.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return The number of checks \code{ng.check}, the number of no checks \code{ng},
#' the number of missing values for checks \code{nmis.check}), the number of
#' missing values for no checks \code{nmis}, the number \code{n.check.0}
#' and list \code{check.0} of checks without data, the number \code{n.check.1}
#' and list \code{check.1} of checks with only one datum, the number of checks
#' with at least two data \code{n.check.2}, the number of replications \code{nr},
#' and the number of rows in the data frame with missing values for factors
#' (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @details This function checks the frequencies and number of missing values in an ABD.
#' For an ANOVA in an ABD it is needed at least two checks with at least 2 valid cases each.
#' @export

ck.abd <- function(trait, geno, rep, data) {
  
  # Check and delete rows with missing values for factors
  
  nmis.fact <- nrow(data[is.na(data[, geno]) | is.na(data[, rep]), ])
  
  if (nmis.fact > 0)
    data <- data[!(is.na(data[, geno]) | is.na(data[, rep])), ]
  
  # Number of replications
  
  nr <- length(unique(data[, rep]))

  # Identify checks and no checks
  
  tfreq <- data.frame(table(data[, geno]))
  check <- as.character(tfreq[tfreq$Freq > 1, 1])
  nocheck <- as.character(tfreq[tfreq$Freq == 1, 1])
  
  # Number of checks and no checks
  
  ng.check <- length(check)
  ng <- length(nocheck)
  
  # Number of missing values for no checks
  
  temp <- data[!(data[, geno] %in% check), ]
  nmis <- sum(is.na(temp[, trait]))
  
  # Number of missing values for checks
  # Factor format to preserve levels in the table of frequencies
  
  temp <- data[data[, geno] %in% check, ]
  temp$geno <- factor(temp$geno)

  temp <- temp[!is.na(temp[, trait]), ]
  tfreq <- table(temp[, geno], temp[, rep])
  nmis.check <- sum(tfreq == 0)

  # Number of checks without data, with 1, and more data

  tfreq <- data.frame(table(temp[, geno]))
  
  n.check.0 <- sum(tfreq$Freq == 0)
  n.check.1 <- sum(tfreq$Freq == 1)
  n.check.2 <- sum(tfreq$Freq > 1)
  
  # List of checks witout data or only one datum
  
  check.0 <- NULL
  if (n.check.0 > 0)
    check.0 <- tfreq[tfreq$Freq == 0, 1]
  
  check.1 <- NULL
  if (n.check.1 > 0)
    check.1 <- tfreq[tfreq$Freq == 1, 1]

  # Return
  
  list(ng.check = ng.check, ng = ng, nmis.check = nmis.check, nmis = nmis, nr = nr,
       n.check.0 = n.check.0, check.0 = check.0, n.check.1 = n.check.1,
       check.1 = check.1, n.check.2 = n.check.2, nmis.fact = nmis.fact)
  
}

#' Check data for a CRD
#'
#' This function checks the frequencies of genotypes in a CRD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param data The name of the data frame.
#' @return Two control values (\code{c1} and \code{c2},
#' the number of genotypes (\code{ng}), the number of replications (\code{nr}),
#' a table with frequencies of valid cases for each genotype, and the number of
#' rows in the data frame with missing values for factors (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication in a CRD and
#' if there is any genotype without data.
#' @export

ck.crd <- function(trait, geno, data) {
  
  # Check and delete rows with missing values for factors
  
  nmis.fact <- nrow(data[is.na(data[, geno]), ])
  
  if (nmis.fact > 0)
    data <- data[!(is.na(data[, geno])), ]

  # Genotypes as factor to preserve levels in the table of frequencies
  
  data[, geno] <- factor(data[, geno])

  # Check frequencies by genotype
  
  subdata <- subset(data, !is.na(data[, trait]))
  tfreq <- table(subdata[, geno])

  # Number of levels
  
  ng <- nlevels(data[, geno])
  nr <- max(tfreq)

  # Controls
  
  c1 <- 1 # Check for zeros. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replication

  if (min(tfreq) == 0)
    c1 <- 0 # State 0: there are zeros
  if (nr > 1)
    c2 <- 1 # State 1: more than one replication

  # Return
  
  list(c1 = c1, c2 = c2, ng = ng, nr = nr, tfreq = tfreq, nmis.fact = nmis.fact)
  
}

#' Check data for a RCBD
#'
#' This function checks the frequencies of genotypes in a RCBD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Four control values (\code{c1}, \code{c2}, \code{c3}, and \code{c4}),
#' the number of missing values \code{nmis}, the proportion of missing values
#' (\code{pmis}), the number of genotypes (\code{ng}), the number of replications
#' (\code{nr}), a table with frequencies of valid cases for each genotype, and
#' the number of rows in the data frame with missing values for factors
#' (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication in a RCBD,
#' if there is any genotype without data or with more data than replications, and
#' if the design is balanced.
#' @export

ck.rcbd <- function(trait, geno, rep, data) {
  
  # Check and delete rows with missing values for factors
  
  nmis.fact <- nrow(data[is.na(data[, geno]) | is.na(data[, rep]), ])
  
  if (nmis.fact > 0)
    data <- data[!(is.na(data[, geno]) | is.na(data[, rep])), ]
  
  # Everything as factor to preserve levels in the table of frequencies
  
  data[, geno] <- factor(data[, geno])
  data[, rep] <- factor(data[, rep])
  
  # Number of levels
  
  ng <- nlevels(data[, geno])
  nr <- nlevels(data[, rep])
  
  # Check frequencies by geno
  
  subdata <- subset(data, !is.na(data[, trait]))
  tfreq <- table(subdata[, geno], subdata[, rep])
  nmis <- sum(tfreq == 0)
  pmis <- mean(tfreq == 0)
  
  # Controls (1 good, 0 bad)
  
  c1 <- 1 # Check for zeros. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replication
  c3 <- 1 # Check for genotypes with more than one datum in a replication
  c4 <- 1 # Check for missing values. Initial state no missing values
  
  if (min(table(subdata[, geno])) == 0)
    c1 <- 0 # State 0: there are zeros
  if (nr > 1)
    c2 <- 1 # State 1: more than one replication
  if (max(tfreq) > 1)
    c3 <- 0 # State 0: some genotypes with addional data
  if (min(tfreq) == 0)
    c4 <- 0 # State 0: missing values
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, nmis = nmis, pmis = pmis,
       ng = ng, nr = nr, tfreq = tfreq, nmis.fact = nmis.fact)
  
}

#' Check data for a full factorial
#'
#' This function checks the frequencies for a full factorial.
#' @param trait The trait to analyze.
#' @param factors The factors.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Four control values (\code{c1}, \code{c2}, \code{c3}, and \code{c4}),
#' the number of missing values \code{nmis}, the proportion of missing values
#' (\code{pmis}), the number of factors (\code{nf}), the number of levels of each
#' factor (\code{nl}), the number of replications (\code{nr}), a table with
#' frequencies of valid cases for each combination of the levels of the factors
#' (\code{tfreq}), a table with frequencies of valid cases for each combination
#' of the levels of the factors in each replication (\code{tfreqr}), and the
#' number of rows in the data frame with missing values for factors (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @details This function checks if there is more than one replication, if there is
#' any combination of the levels of the factors without data or with more data than
#' replications, and if the design is balanced.
#' @export

ck.f <- function(trait, factors, rep, data) {
  
  # Check and delete rows with missing values for factors
  
  cond <- apply(data[, c(factors, rep)], 1, function(x) sum(is.na(x)) > 0)
  
  nmis.fact <- sum(cond)
  
  if (nmis.fact > 0)
    data <- data[!cond, ]
  
  # Number of factors
  
  nf <- length(factors)
  
  # Number of levels for factors
  
  nl <- apply(data[, factors], 2, function(x) length(unique(x)))
  
  # Number of replications

  nr <- length(unique(data[, rep]))

  # Check frequencies
  
  nmis <- sum(is.na(data[, trait]))
  pmis <- mean(is.na(data[, trait]))
  
  subdata <- subset(data, !is.na(data[, trait]))
  
  expr <- 'table(subdata[, factors[1]]'
  
  for (i in 2:nf)
    expr <- paste0(expr, ', subdata[, factors[', i, ']]')
  
  expr1 <- paste0(expr, ')')
  expr2 <- paste0(expr, ', subdata[, rep])')
  
  tfreq <- eval(parse(text = expr1))
  tfreqr <- eval(parse(text = expr2))

  # Controls
  
  c1 <- 1 # Check for zero frequencies. Initial state no zeros
  c2 <- 0 # Check for replicates. Initial state only one replicate
  c3 <- 1 # Check for genotypes with more than one datum in a replication of one environment
  c4 <- 1 # Check for missing values. Initial state no missing values
    
  if (min(tfreq) == 0)
    c1 <- 0 # State 0: there are zeros
  if (nr > 1)
    c2 <- 1 # State 1: more than one replicate
  if (max(tfreqr) > 1)
    c3 <- 0 # State 0: some genotypes with addional data
  if (min(tfreqr) == 0)
    c4 <- 0 # State 0: missing values
    
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, nmis = nmis, pmis = pmis, nf = nf,
       nl = nl, nr = nr, tfreq = tfreq, tfreqr = tfreqr, nmis.fact = nmis.fact)
  
}

#' Check data for a Wescott layout
#'
#' This function checks the grid of checks on the Wescott layout and
#' the number of missing values.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param data The name of the data frame.
#' @return Four control values (\code{c1}, \code{c2}, \code{c3}, and \code{c4},
#' for the grid of checks, the number of missing values for checks (\code{nmis.check})
#' and genotypes \code{nmis}, the proportion of missing values for checks
#' (\code{pmis.check}) and genotypes (\code{pmis}), and the number of rows
#' in the data frame with missing values for factors (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @details This function checks the grid of checks for the Wescoot layout and
#' calculates the number of missing values.
#' @export

ck.w <- function(trait, geno, ch1, ch2, row, col, ncb, data) {
  
  # Check and delete rows with missing values for factors
  
  nmis.fact <- nrow(data[is.na(data[, geno]) | is.na(data[, row]) | is.na(data[, col]), ])
  
  if (nmis.fact > 0)
    data <- data[!(is.na(data[, geno]) | is.na(data[, rep]) | is.na(data[, col])), ]
  # Checks
  
  checks <- c(ch1, ch2)
  
  # Rows and columns as numbers
  
  data[, row] <- as.numeric(data[, row])
  data[, col] <- as.numeric(data[, col])

  # Number of rows and columns
  
  nr.min <- min(data[, row])
  nc.min <- min(data[, col])
  nc.max <- max(data[, col])
  
  # Controls
  
  c1 <- 0 # All column checks with checks
  c2 <- 0 # All column genotypes with genotypes
  c3 <- 0 # Alternating checks without in correlative row order
  c4 <- 0 # All genotypes with checks to the left and right

  # Columns with checks
  
  cch <- seq(nc.min, nc.max, ncb + 1)
  
  # Check columns with checks
  
  if (sum(!(data[data[, col] %in% cch, geno] %in% checks)) > 0)
    c1 <- 1
  
  # Check columns with genotypes
  
  if (sum(data[!(data[, col] %in% cch), geno] %in% checks) > 0)
    c2 <- 1
  
  # Alternating checks
  
  for (i in cch)
    for (j in (min(data[data[, col] == i, row]) + 1):max(data[data[, col] == i, row]))
      if (data[data[, col] == i & data[, row] == j, geno] == data[data[, col] == i & data[, row] == j - 1, geno])
        c3 <- 1
  
  for (i in 2:length(cch))
    if (data[data[, col] == cch[i] & data[, row] == nr.min, geno] == data[data[, col] == cch[i - 1] & data[, row] == nr.min, geno])
      c3 <- 1
  
  # All genotypes must have one check to the left and one to the right
  
  for (i in 1:dim(data)[1]) {
    rows <- data[i, row]
    columns <- (data[i, col] - ncb):(data[i, col] + ncb)
    temp <- data[data[, row] == rows & data[, col] %in% columns, geno]
    if (sum(temp %in% checks) == 0)
      c4 <- 1
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
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, nmis = nmis, pmis = pmis,
       nmis.check = nmis.check, pmis.check = pmis.check, nmis.fact = nmis.fact)
  
}
