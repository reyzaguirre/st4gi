#' Check number of genotypes and replications
#' 
#' This function cheks the number of genotypes and replications for different designs.
#' 
#' @param geno The genotypes
#' @param design The statistical design
#' @param dfr The name of the data frame
#' @return The number of genotypes (\code{ng}), the number of checks (\code{ng.check}),
#' the list of genotypes (\code{lg}), the list of checks (\code{lg.check}),
#' and the number of replications (\code{nr}).
#' @author Raul Eyzaguirre.
#' @examples 
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check the design
#' ck.gr('geno', 'block', 'rcbd', dfr = dfr)
#' @export

ck.gr <- function(geno, rep = NULL, design = c('crd', 'rcbd', 'abd'), dfr) {

  # match arguments
  
  design <- match.arg(design)
  
  # Check and remove rows with missing values for factors
  
  dfr <- rm.fna(c(geno, rep), dfr)$dfr

  # Number of genotypes and checks
  
  if (design == 'abd') {
    tfreq <- data.frame(table(dfr[, geno]))
    lg.check <- as.character(tfreq[tfreq$Freq > 1, 1])
    lg <- as.character(tfreq[tfreq$Freq == 1, 1])
    ng.check <- length(lg.check)
    ng <- length(lg)
  } else {
    lg <- as.character(unique(dfr[, geno]))
    ng <- length(lg)
    lg.check <- NULL
    ng.check <- NULL
  }
  
  # Number of replications
  
  if (design %in% c('abd', 'rcbd')) {
    nr <- length(unique(dfr[, rep]))
  } else {
    tfreq <- table(dfr[, geno])
    nr <- max(tfreq)
  }
  
  # Return
  
  list(ng = ng, ng.check = ng.check, lg = lg, lg.check = lg.check, nr = nr)
  
}

#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param dfr The name of the data frame.
#' @return For each replication, the number (\code{nplot}) and list (\code{lplot})
#' of plots (unique row and column position) with more than one genotype, and the 
#' number (\code{nr}) and list (\code{lr}) of replications.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check positions
#' ck.pos('row', 'col', 'block', dfr = dfr)
#' @export

ck.pos <- function(row, col, rep = NULL, dfr) {
  
  # Define replications
  
  if (is.null(rep)) {
    dfr[, "rep"] <- 1
    rep <- "rep"
  }
  
  # Check and remove rows with missing values for factors
  
  dfr <- rm.fna(c(row, col, rep), dfr)$dfr
  
  # Number of replications and levels

  lr <- sort(unique(dfr[, rep]))
  nr <- length(lr)
  
  # Number and list of plots with more than one genotype
  
  nplot <- NULL
  lplot <- list()
  
  # Check row and column

  for (i in 1:nr) {
    
    # Compute frequencies
    
    temp <- dfr[dfr[, rep] == lr[i], ]
    ttt <- as.data.frame(table(temp[, row], temp[, col]))
    colnames(ttt) <- c('Row', 'Column', 'Freq')
    
    # Number of plots with problems
    
    nplot[i] <- dim(ttt[ttt$Freq > 1, ])[1]
    
    # List of plots with problems if any
    
    if (nplot[i] > 0)
      lplot[[i]] <- ttt[ttt$Freq > 1, ] 
  }
  
  # Return
  
  list(nplot = nplot, lplot = lplot, nr = nr, lr = lr)
  
}

#' Check data for an ABD
#'
#' This function checks the frequencies of genotypes in an ABD.
#' @param trait The trait to analyze.
#' @param geno The genotypes including checks.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return The number of checks \code{ng.check}, the number of no checks \code{ng},
#' the number of missing values for checks \code{nmis.check}), the number of
#' missing values for no checks \code{nmis}, the number \code{ncheck.0}
#' and list \code{check.0} of checks without data, the number \code{ncheck.1}
#' and list \code{check.1} of checks with only one datum, the number of checks
#' with at least two data \code{ncheck.2}, the number of replications \code{nr},
#' and the number of rows in the data frame with missing values for factors
#' (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples
#' # Create design
#' dfr <- cr.abd(1:50, c('a', 'b', 'd'), 5, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(65)
#' # Delete some values
#' dfr[c(1, 5, 7, 56), 'y'] <- NA
#' # Delete some values for classification factors
#' dfr[64, 'geno'] <- NA
#' # Check the design
#' ck.abd('y', 'geno', 'block', dfr)
#' @export

ck.abd <- function(trait, geno, rep, dfr) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(c(geno, rep), dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac

  # Number and list of checks and nochecks, and number of replications
  
  out <- ck.gr(geno, rep, 'abd', dfr)
  ng <- out$ng
  lg <- out$lg
  ng.check <- out$ng.check
  lg.check <- out$lg.check
  nr <- out$nr
  
  # Number of missing values for no checks
  
  temp <- dfr[dfr[, geno] %in% lg, ]
  nmis <- sum(is.na(temp[, trait]))
  
  # Number of missing values for checks
  # Factor format to preserve levels in the table of frequencies
  
  temp <- dfr[dfr[, geno] %in% lg.check, ]
  temp[, geno] <- factor(temp[, geno])
  temp[, rep] <- factor(temp[, rep])
  temp <- temp[!is.na(temp[, trait]), ]
  tfreq <- table(temp[, geno], temp[, rep])
  nmis.check <- sum(tfreq == 0)

  # Number of checks without data, with 1, and more data

  tfreq <- data.frame(table(temp[, geno]))
  
  ncheck.0 <- sum(tfreq$Freq == 0)
  ncheck.1 <- sum(tfreq$Freq == 1)
  ncheck.2 <- sum(tfreq$Freq > 1)
  
  # List of checks witout data or only one datum
  
  check.0 <- NULL
  if (ncheck.0 > 0)
    check.0 <- tfreq[tfreq$Freq == 0, 1]
  
  check.1 <- NULL
  if (ncheck.1 > 0)
    check.1 <- tfreq[tfreq$Freq == 1, 1]

  # Return
  
  list(ng.check = ng.check, ng = ng, nmis.check = nmis.check, nmis = nmis,
       ncheck.0 = ncheck.0, check.0 = check.0, ncheck.1 = ncheck.1,
       check.1 = check.1, ncheck.2 = ncheck.2, nmis.fac = nmis.fac, nr = nr)
  
}

#' Check data for a CRD
#'
#' This function checks the frequencies of genotypes in a CRD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param dfr The name of the data frame.
#' @return The number of genotypes (\code{ng}), the number of genotypes without
#' data (\code{ng.0}), the number of replications (\code{nr}), a table with
#' frequencies of valid cases for each genotype (\code{tfreq}), and the number
#' of rows in the data frame with missing values for factors (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples
#' # Create design
#' dfr <- cr.crd(1:50, 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(150)
#' # Delete some values
#' dfr[c(1, 5, 56, 77, 111), 'y'] <- NA
#' # Delete some values for classification factors
#' dfr[c(27, 48), 'geno'] <- NA
#' # Check the design
#' ck.crd('y', 'geno', dfr)
#' @export

ck.crd <- function(trait, geno, dfr) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(geno, dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac

  # Number of genotypes and replications
  
  out <- ck.gr(geno, NULL, 'crd', dfr)
  ng <- out$ng
  nr <- out$nr

  # Genotypes as factor to preserve levels in the table of frequencies
  
  temp <- dfr
  temp[, geno] <- factor(temp[, geno])
  
  # Frequencies for genotypes
  
  temp <- temp[!is.na(temp[, trait]), ]
  tfreq <- table(temp[, geno])

  # Number of genotypes without data
  
  ng.0 <- sum(tfreq == 0)
  
  # Return
  
  list(ng.0 = ng.0, ng = ng, nr = nr, tfreq = tfreq, nmis.fac = nmis.fac)
  
}

#' Check data for a RCBD
#'
#' This function checks the frequencies of genotypes in a RCBD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return The number of genotypes without data (\code{ng.0}), the number of
#' genotypes with more than one plot in a given block (\code{ng.2}), the number
#' of missing values \code{nmis}, the proportion of missing values (\code{pmis}),
#' the number of genotypes (\code{ng}), the number of replications (\code{nr}),
#' a table with frequencies of valid cases for each genotype, and the number of
#' rows in the data frame with missing values for factors (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples 
#' # Create design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(60)
#' # Delete some values
#' dfr[c(1, 5, 16, 17), 'y'] <- NA
#' # Check the design
#' ck.rcbd('y', 'geno', 'block', dfr)
#' @export

ck.rcbd <- function(trait, geno, rep, data) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(c(geno, rep), dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac
  
  # Number of genotypes and replications
  
  out <- ck.gr(geno, rep, 'rcbd', dfr)
  ng <- out$ng
  nr <- out$nr

  # Genotypes and replications as factors to preserve levels in the table of frequencies
  
  temp <- dfr
  temp[, geno] <- factor(temp[, geno])
  temp[, rep] <- factor(temp[, rep])

  # Frequencies by geno
  
  temp <- temp[!is.na(temp[, trait]), ]
  tfreq <- table(temp[, geno], temp[, rep])
  
  # Number of missing values
  
  nmis <- sum(tfreq == 0)
  pmis <- mean(tfreq == 0)
  
  # Number of genotypes with more than one plot in a given block
  
  ng.2 <- sum(tfreq > 1)
  
  # Number of genotypes without data
  
  tfreq <- table(temp[, geno])
  ng.0 <- sum(tfreq == 0)

  # Return
  
  list(ng.0 = ng.0, ng.2 = ng.2, nmis = nmis, pmis = pmis,
       ng = ng, nr = nr, tfreq = tfreq, nmis.fac = nmis.fac)
  
}

#' Check data for a full factorial
#'
#' This function checks the frequencies for a full factorial.
#' @param trait The trait to analyze.
#' @param factors The factors.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return control values (\code{c1}, and \code{c3}, the number
#' of missing values \code{nmis}, the proportion of missing values (\code{pmis}),
#' the number of factors (\code{nf}), the number of levels of each factor (\code{nl}),
#' the number of replications (\code{nr}), a table with frequencies of valid cases
#' for each combination of the levels of the factors (\code{tfreq}), a table with
#' frequencies of valid cases for each combination of the levels of the factors in
#' each replication (\code{tfreqr}), and the number of rows in the data frame with
#' missing values for factors (\code{nmis.fact}).
#' @author Raul Eyzaguirre.
#' @export

ck.f <- function(trait, factors, rep, dfr) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(c(factors, rep), dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac
  
  # Number of factors
  
  nf <- length(factors)
  
  # Number of levels for factors
  
  nl <- apply(dfr[, factors], 2, function(x) length(unique(x)))
  
  # Number of replications

  nr <- length(unique(dfr[, rep]))

  # Number of missing values
  
  nmis <- sum(is.na(dfr[, trait]))
  pmis <- mean(is.na(dfr[, trait]))
  
  # Factors and replications as factors to preserve levels in the table of frequencies
  
  temp <- dfr
  for (i in 1:nf)
    temp[, factors[i]] <- factor(temp[, factors[i]])
  temp[, rep] <- factor(temp[, rep])

  # Calculate frequencies
  
  temp <- temp[!is.na(temp[, trait]), ]
  
  expr <- 'table(temp[, factors[1]]'
  
  for (i in 2:nf)
    expr <- paste0(expr, ', temp[, factors[', i, ']]')
  
  expr1 <- paste0(expr, ')')
  expr2 <- paste0(expr, ', temp[, rep])')
  
  tfreq <- eval(parse(text = expr1))
  tfreqr <- eval(parse(text = expr2))

  # Controls

  c1 <- 1 # State 1: No zeros
  c3 <- 1 # State 1: Each genotype only once in each replication in each environment

  if (min(tfreq) == 0)
    c1 <- 0 # State 0: There are zeros
  if (max(tfreqr) > 1)
    c3 <- 0 # State 0: Some genotypes with addional data

  # Return
  
  list(c1 = c1, c3 = c3, nmis = nmis, pmis = pmis, nf = nf, nl = nl,
       nr = nr, tfreq = tfreq, tfreqr = tfreqr, nmis.fact = nmis.fact)
  
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
