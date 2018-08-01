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
