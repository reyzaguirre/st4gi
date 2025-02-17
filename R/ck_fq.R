#' Check frequencies
#' 
#' This function cheks the frequencies of valid cases for treatments and replications.
#' @param dfr The name of the data frame.
#' @param y The name of the column for the variable to analyze.
#' @param factors The names of the columns that identify the factors.
#' @param rep The name of the column that identifies the replications,
#' default is \code{NULL} for a CRD.
#' @return A table of frequencies of valid cases for all factors' levels combinations
#' (\code{tf}), a table of frequencies of valid cases for all factors' levels and
#' replications combinations (\code{tfr}), the number of missing values \code{nmis},
#' and the proportion of missing values (\code{pmis}).
#' @author Raul Eyzaguirre.
#' @examples 
#' ## Example 1
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(60)
#' # Delete some values
#' dfr[c(1, 5, 16, 17), 'y'] <- NA
#' # Check the frequencies
#' ck.fq(dfr, "y", "geno", "block")
#' 
#' ## Example 2
#' # Create a design
#' A <- paste0("a", 1:5)
#' B <- paste0("b", 1:3)
#' dfr <- cr.f(c("A", "B"), list(A, B), "rcbd", 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(45)
#' # Delete some values
#' dfr[c(5, 10, 24), 'y'] <- NA
#' # Check the frequencies
#' ck.fq(dfr, "y", c("A", "B"), "block")
#' @export

ck.fq <- function(dfr, y, factors, rep = NULL) {
  
  # Number of missing values
  
  nmis <- sum(is.na(dfr[, y]))
  pmis <- mean(is.na(dfr[, y]))
  
  # tfr NULL for rep = NULL
  
  tfr <- NULL
  
  # Number of factors
  
  nf <- length(factors)
  
  # Factors and replications as factors to preserve levels in the table of frequencies
  
  tmp <- dfr
  for (i in 1:nf)
    tmp[, factors[i]] <- factor(tmp[, factors[i]])
  if (!is.null(rep))
    tmp[, rep] <- factor(tmp[, rep])
  
  # Calculate frequencies
  
  tmp <- tmp[!is.na(tmp[, y]), ]

  if (nf == 1) {
    tf <- table(tmp[, factors])
    if (!is.null(rep))
      tfr <- table(tmp[, factors], tmp[, rep])
  } else {
    expr <- 'table(tmp[, factors[1]]'
    for (i in 2:nf)
      expr <- paste0(expr, ', tmp[, factors[', i, ']]')
    expr1 <- paste0(expr, ')')
    tf <- eval(parse(text = expr1))
    if (!is.null(rep)) {
      expr2 <- paste0(expr, ', tmp[, rep])')
      tfr <- eval(parse(text = expr2))
    }
  }
  
  # Return
  
  list(tf = tf, tfr = tfr, nmis = nmis, pmis = pmis)
  
}
