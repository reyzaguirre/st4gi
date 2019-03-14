#' Check factors structure
#' 
#' This function cheks the structure of factors.
#' @param factors The factors.
#' @param rep The replications, \code{NULL} for a CRD.
#' @param dfr The name of the data frame.
#' @return The number of factors (\code{nf}), the number of levels of the factors
#' (\code{nl}), the lists of levels of factors (\code{lf}), the number of treatments
#' (\code{nt}), the list of treatments (\code{lt}), the number of replications
#' (\code{nrep}), the list of replications (\code{lrep}), the number of rows
#' with missing values for factors (\code{nmis.fac}), and the data frame after
#' removal of all these rows.
#' @author Raul Eyzaguirre.
#' @examples 
#' ## Example 1
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check the design
#' ck.fs("geno", "block", dfr)
#' 
#' ## Example 2
#' # Create a design
#' A <- paste0("a", 1:5)
#' B <- paste0("b", 1:3)
#' dfr <- cr.f(c("A", "B"), list(A, B), "rcbd", 3, 10)
#' dfr <- dfr$book
#' # Check the design
#' ck.fs(c("A", "B"), "block", dfr)
#' @export

ck.fs <- function(factors, rep, dfr) {

  # Check missing values for factors
  
  if (is.null(rep)) {
    cond <- apply(data.frame(dfr[, factors]), 1, function(x) sum(is.na(x)) > 0)
  } else {
    cond <- apply(dfr[, c(factors, rep)], 1, function(x) sum(is.na(x)) > 0)
  }
  
  # Number of missing values for factors
  
  nmis.fac <- sum(cond)
  
  # Remove rows with missing values for factors
  
  if (nmis.fac > 0)
    dfr <- dfr[!cond, ]

  # Number of factors
  
  nf <- length(factors)
  
  # Number of levels for factors
  
  nl <- apply(data.frame(dfr[, factors]), 2, function(x) length(unique(x)))
  names(nl) <- factors
  
  # Levels of factors
  
  lf <- apply(data.frame(dfr[, factors]), 2, function(x) unique(x))
  if (nf == 1 & length(lf) > 1)
    colnames(lf) <- factors
  if (nf == 1 & length(lf) == 1)
    names(lf) <- factors

  # Define treatments
  
  treat <- factors[1]
  if (nf > 1)
    for (i in 2:nf)
      treat <- paste0(treat, factors[i])

  dfr[, treat] <- dfr[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf)
      dfr[, treat] <- paste(dfr[, treat], dfr[, factors[i]], sep = "_")

  # Number of treatments
  
  lt <- sort(unique(dfr[, treat]))
  nt <- length(lt)

  # Number and levels of replications
  
  if (is.null(rep)) {
    tfreq <- table(dfr[, treat])
    lrep <- NULL
    nrep <- max(tfreq)
  } else {
    lrep <- sort(unique(dfr[, rep]))
    nrep <- length(lrep)
  }
  
  # Return
  
  list(nf = nf, nl = nl, lf = lf, nt = nt, lt = lt, nrep = nrep, lrep = lrep,
       nmis.fac = nmis.fac, dfr = dfr)
  
}
