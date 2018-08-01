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
#' @examples 
#' Create a design
#' A <- paste0("a", 1:5)
#' B <- paste0("b", 1:3)
#' dfr <- cr.f(c("A", "B"), list(A, B), "rcbd", 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(45)
#' # Delete some values
#' dfr[c(4, 5, 12), 'y'] <- NA
#' # Check the design
#' ck.f("y", c("A", "B"), "block", dfr)
#' @export

ck.f <- function(trait, factors, rep, dfr) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(c(factors, rep), dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac
  
  # Number of factors, levels and replications
  
  out <- ck.fs(factors, rep, 'rcbd', dfr)
  nf <- out$nf
  nl <- out$nl
  nr <- out$nr

  # Number of missing values
  
  nmis <- sum(is.na(dfr[, trait]))
  pmis <- mean(is.na(dfr[, trait]))

  # Frequencies for factors and replications
  
  out <- ck.fq(trait, factors, rep, dfr)
  tf <- out$tf
  tfr <- out$tfr

  # Number of treatments without data
  
  nt.0 <- sum(tf == 0)
  
  # Number of treatments that appear more than once in a given replication
  
  nt.mult <- sum(tfr > 1)

  # Return
  
  list(nt.0 = nt.0, nt.mult = nt.mult, nmis = nmis, pmis = pmis, nf = nf, nl = nl,
       nr = nr, nmis.fac = nmis.fac)
  
}
