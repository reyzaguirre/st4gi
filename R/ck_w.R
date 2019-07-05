#' Check data for a Wescott layout
#'
#' This function checks the grid of checks on the Wescott layout and
#' the number of missing values.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param ck1 Name of check 1.
#' @param ck2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param dfr The name of the data frame.
#' @return Four control values (\code{c1}, \code{c2}, \code{c3}, and \code{c4},
#' for the grid of checks, the number of missing values for checks (\code{nmis.ck})
#' and genotypes \code{nmis}, the proportion of missing values for checks
#' (\code{pmis.ck}) and genotypes (\code{pmis}), and the number of rows
#' in the data frame with missing values for factors (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.w(1:1000, "A", "B", 50, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(1125)
#' # Delete some values
#' dfr[c(11, 165, 569, 914), 'y'] <- NA
#' # Check the design
#' ck.w("y", "geno", "A", "B", "row", "col", 10, dfr)
#' @export

ck.w <- function(trait, geno, ck1, ck2, row, col, ncb, dfr) {
  
  # Check and remove rows with missing values for factors
  
  out <- ck.fs(c(geno, row, col), NULL, dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac

  # Checks
  
  checks <- c(ck1, ck2)
  
  # Rows and columns as numbers
  
  dfr[, row] <- as.numeric(dfr[, row])
  dfr[, col] <- as.numeric(dfr[, col])

  # Number of rows and columns
  
  nr.min <- min(dfr[, row])
  nc.min <- min(dfr[, col])
  nc.max <- max(dfr[, col])
  
  # Controls
  
  c1 <- 0 # All column checks with checks (no genotypes)
  c2 <- 0 # All column genotypes with genotypes (no checks)
  c3 <- 0 # All column checks with checks (no missing plots)
  c4 <- 0 # Alternating checks in rows and columns
  c5 <- 0 # All genotypes with checks to the left and right

  # Columns with checks
  
  cck <- seq(nc.min, nc.max, ncb + 1)
  
  # Check columns with checks
  
  if (sum(!(dfr[dfr[, col] %in% cck, geno] %in% checks)) > 0)
    c1 <- 1
  
  # Check columns with genotypes
  
  if (sum(dfr[!(dfr[, col] %in% cck), geno] %in% checks) > 0)
    c2 <- 1
  
  # Alternating checks on columns
  
  for (i in cck)
    for (j in (min(dfr[dfr[, col] == i, row]) + 1):max(dfr[dfr[, col] == i, row])) {
      val1 <- dfr[dfr[, col] == i & dfr[, row] == j, geno]
      val2 <- dfr[dfr[, col] == i & dfr[, row] == j - 1, geno]
      if (length(val1) == 0 | length(val2) == 0) {
        c3 <- 1
      } else {
        if (val1 == val2)
          c4 <- 1
      }
    }
  
  # Alternating checks on rows

  for (i in 2:length(cck)) {
    val1 <- dfr[dfr[, col] == cck[i] & dfr[, row] == nr.min, geno]
    val2 <- dfr[dfr[, col] == cck[i - 1] & dfr[, row] == nr.min, geno]
    if (length(val1) == 1 & length(val2) == 1)
      if (val1 == val2)
        c4 <- 1
  }
  
  # All genotypes must have one check to the left and one to the right
  
  for (i in 1:dim(dfr)[1]) {
    rows <- dfr[i, row]
    columns <- (dfr[i, col] - ncb):(dfr[i, col] + ncb)
    temp <- dfr[dfr[, row] == rows & dfr[, col] %in% columns, geno]
    if (sum(temp %in% checks) == 0)
      c5 <- 1
  }

  # Number of missing values for checks
  
  temp <- dfr[dfr[, col] %in% cck, ]
  nmis.ck <- sum(is.na(temp[, trait]))
  pmis.ck <- mean(is.na(temp[, trait]))

  # Number of missing values for genotypes
  
  temp <- dfr[!(dfr[, col] %in% cck), ]
  nmis <- sum(is.na(temp[, trait]))
  pmis <- mean(is.na(temp[, trait]))
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, nmis = nmis, pmis = pmis,
       nmis.ck = nmis.ck, pmis.ck = pmis.ck, nmis.fac = nmis.fac)
  
}
