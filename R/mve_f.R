#' Estimation of missing values for a factorial experiment
#'
#' Function to estimate missing values for factorial experiment with a CRD
#' or a RCBD by the least squares method.
#' @param trait The trait to estimate missing values.
#' @param factors The factors.
#' @param rep The replications or blocks, \code{NULL} for a CRD.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated
#' values.
#' @author Raul Eyzaguirre.
#' @examples
#' # A data frame with some missing values
#' temp <- asc
#' temp$dm[c(3, 10, 115)] <- NA
#' # Estimate the missing values
#' mve.f("dm", c("geno", "treat"), NULL, temp)
#' @export

mve.f <- function(trait, factors, rep, dfr, maxp = 0.1, tol = 1e-06) {
  
  # Check data
  
  lc <- ck.f(trait, factors, rep, dfr)
  
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
  
  # Get a copy of trait for estimated values and for temporary values
  
  trait.est <- paste0(trait, ".est")
  dfr[, trait.est] <- dfr[, trait]
  dfr[, "ytemp"] <- dfr[, trait]
  
  # Create expression for list of factors
  
  lf.expr <- 'list(dfr[, factors[1]]'
  
  for (i in 2:lc$nf)
    lf.expr <- paste0(lf.expr, ', dfr[, factors[', i, ']]')
  
  lf.expr <- paste0(lf.expr, ')')
  
  # Compute means over replications

  tmeans <- tapply(dfr[, trait], eval(parse(text = lf.expr)), mean, na.rm = TRUE)

  # Store means in ytemp for missing values
  
  for (i in 1:length(dfr[, trait]))
    if (is.na(dfr[i, trait])) {
      expr <- paste0('tmeans[dfr[', i, ', factors[1]]')
      for (j in 2:lc$nf)
        expr <- paste0(expr, ', dfr[', i, ', factors[', j, ']]')
      expr <- paste0(expr, ']')
      dfr[i, "ytemp"] <- eval(parse(text = expr))
    }
  
  # Estimate missing values for a crd
  
  if (is.null(rep))
    dfr[, trait.est] <- dfr[, "ytemp"]
  
  # Estimate missing values for a rcbd
  
  if (!is.null(rep)) {
    
    # Total number of treatments
    
    nt <- prod(lc$nl)
    
    # Vectors of estimated missing values to check convergence
    
    emv1 <- array(0, lc$nmis)
    emv2 <- array(0, lc$nmis)
    
    # Value to control convergence
    
    cc <- max(dfr[, trait], na.rm = TRUE)
    
    # Maximum number of iteration control
    
    cont <- 0
    
    # Iterate
    
    while (cc > max(dfr[, trait], na.rm = TRUE) * tol & cont < 100) {
      
      cont <- cont + 1
      
      for (i in 1:length(dfr[, trait]))
        if (is.na(dfr[i, trait])) {
      
          # Compute sums
          
          dfr[i, "ytemp"] <- NA
          sum1 <- tapply(dfr[, "ytemp"], eval(parse(text = lf.expr)), sum, na.rm = TRUE)
          sum2 <- tapply(dfr[, "ytemp"], dfr[, rep], sum, na.rm = TRUE)
          sum3 <- sum(dfr[, "ytemp"], na.rm = TRUE)
          
          # Get estimate
          
          expr <- paste0('sum1[dfr[', i, ', factors[1]]')
          for (j in 2:lc$nf)
            expr <- paste0(expr, ', dfr[', i, ', factors[', j, ']]')
          expr <- paste0(expr, ']')
          
          mv.num <- nt * eval(parse(text = expr)) + lc$nrep * sum2[dfr[i, rep]] - sum3
          mv.den <- nt * lc$nrep - nt - lc$nrep + 1

          dfr[i, trait.est] <- mv.num / mv.den
          
          dfr[i, "ytemp"] <- dfr[i, trait.est]
        }
      
      emv1 <- emv2
      emv2 <- dfr[is.na(dfr[, trait]), trait.est]
      cc <- max(abs(emv1 - emv2))
    }
  }
  
  # Return
  
  dfr[, c(factors, rep, trait, trait.est)]
  
}
