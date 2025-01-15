#' Estimation of missing values for a factorial experiment
#'
#' Function to estimate missing values for factorial experiment with a CRD
#' or a RCBD by the least squares method.
#' @param dfr The name of the data frame.
#' @param y The name of the column for the variable to estimate missing values.
#' @param factors The names of the columns that identify the factors.
#' @param rep The name of the column that identifies the replications or blocks,
#' default is \code{NULL} for a CRD.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{y}
#' and \code{y.est} with the original data and the original data plus the estimated
#' values.
#' @author Raul Eyzaguirre.
#' @examples
#' # A data frame with some missing values
#' tmp <- asc
#' tmp$dm[c(3, 10, 115)] <- NA
#' # Estimate the missing values
#' mve.f(tmp, "dm", c("geno", "treat"))
#' @export

mve.f <- function(dfr, y, factors, rep = NULL, maxp = 0.1, tol = 1e-06) {
  
  # Check data
  
  lc <- ck.f(dfr, y, factors, rep)
  
  # Error messages
  
  if (lc$nmis.fac > 0)
    stop("There are missing values for classification factors.")

  if (lc$nt.0 > 0)
    stop("Some factor levels' combinations have zero frequency.")
  
  if (lc$nrep == 1)
    stop("There is only one replication. Inference is not possible with one replication.")
  
  if (lc$nt.mult > 0)
    stop("Some factor levels' combinations have additional replications.")

  if (lc$nmis == 0)
    stop("The data set is balanced. There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste0("Too many missing values (", format(lc$pmis * 100, digits = 3), "%)."))
  
  if (sum(lc$nl < 2) > 0)
    stop("There are factors with only one level. This is is not a factorial experiment.")
  
  # Get a copy of y for estimated values and for temporary values
  
  y.est <- paste0(y, ".est")
  dfr[, y.est] <- dfr[, y]
  dfr[, "ytmp"] <- dfr[, y]
  
  # Create expression for list of factors
  
  lf.expr <- 'list(dfr[, factors[1]]'
  
  for (i in 2:lc$nf)
    lf.expr <- paste0(lf.expr, ', dfr[, factors[', i, ']]')
  
  lf.expr <- paste0(lf.expr, ')')
  
  # Compute means over replications

  tmeans <- tapply(dfr[, y], eval(parse(text = lf.expr)), mean, na.rm = TRUE)

  # Store means in ytmp for missing values
  
  for (i in 1:length(dfr[, y]))
    if (is.na(dfr[i, y])) {
      expr <- paste0('tmeans[dfr[', i, ', factors[1]]')
      for (j in 2:lc$nf)
        expr <- paste0(expr, ', dfr[', i, ', factors[', j, ']]')
      expr <- paste0(expr, ']')
      dfr[i, "ytmp"] <- eval(parse(text = expr))
    }
  
  # Estimate missing values for a crd
  
  if (is.null(rep))
    dfr[, y.est] <- dfr[, "ytmp"]
  
  # Estimate missing values for a rcbd
  
  if (!is.null(rep)) {
    
    # Total number of treatments
    
    nt <- prod(lc$nl)
    
    # Vectors of estimated missing values to check convergence
    
    emv1 <- array(0, lc$nmis)
    emv2 <- array(0, lc$nmis)
    
    # Value to control convergence
    
    cc <- max(dfr[, y], na.rm = TRUE)
    
    # Maximum number of iteration control
    
    cont <- 0
    
    # Iterate
    
    while (cc > max(dfr[, y], na.rm = TRUE) * tol & cont < 100) {
      
      cont <- cont + 1
      
      for (i in 1:length(dfr[, y]))
        if (is.na(dfr[i, y])) {
      
          # Compute sums
          
          dfr[i, "ytmp"] <- NA
          sum1 <- tapply(dfr[, "ytmp"], eval(parse(text = lf.expr)), sum, na.rm = TRUE)
          sum2 <- tapply(dfr[, "ytmp"], dfr[, rep], sum, na.rm = TRUE)
          sum3 <- sum(dfr[, "ytmp"], na.rm = TRUE)
          
          # Get estimate
          
          expr <- paste0('sum1[dfr[', i, ', factors[1]]')
          for (j in 2:lc$nf)
            expr <- paste0(expr, ', dfr[', i, ', factors[', j, ']]')
          expr <- paste0(expr, ']')
          
          mv.num <- nt * eval(parse(text = expr)) + lc$nrep * sum2[dfr[i, rep]] - sum3
          mv.den <- nt * lc$nrep - nt - lc$nrep + 1

          dfr[i, y.est] <- mv.num / mv.den
          
          dfr[i, "ytmp"] <- dfr[i, y.est]
        }
      
      emv1 <- emv2
      emv2 <- dfr[is.na(dfr[, y]), y.est]
      cc <- max(abs(emv1 - emv2))
    }
  }
  
  # Return
  
  dfr[, c(factors, rep, y, y.est)]
  
}
