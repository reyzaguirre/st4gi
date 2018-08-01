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
