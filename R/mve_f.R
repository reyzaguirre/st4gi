#' Estimation of missing values for a factorial experiment
#'
#' Function to estimate missing values for factorial experiment with a CRD
#' or a RCBD by the least squares method.
#' @param trait The trait to estimate missing values.
#' @param factors The factors.
#' @param rep The replications or blocks.
#' @param design The statistical design, \code{crd} or \code{rcbd}.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param tol Tolerance for the convergence of the iterative estimation process.
#' @return It returns a data frame with the experimental layout and columns \code{trait}
#' and \code{trait.est} with the original data and the original data plus the estimated values.
#' @author Raul Eyzaguirre.
#' @details A \code{data.frame} with data for a factorial experiment with at least two
#' replications and at least one datum for each factor levels' combination must be loaded.
#' Experimental data with only one replication, any factor levels' combination without data,
#' or more missing values than specified in \code{maxp} will generate an error message.
#' @examples
#' # The data
#' str(asc)
#' 
#' # A copy with some missing values
#' temp <- asc
#' temp$dm[c(3, 10, 115)] <- NA
#'
#' # Estimate the missing values
#' mve.f("dm", c("geno", "treat"), "rep", "crd", temp)
#' @export

mve.f <- function(trait, factors, rep, design = c("crd", "rcbd"),
                  data, maxp = 0.1, tol = 1e-06) {
  
  # match arguments
  
  design <- match.arg(design)
  
  # Check data
  
  lc <- ck.f(trait, factors, rep, data)
  
  # Error messages
  
  if (lc$nmis.fact > 0)
    stop("There are missing values for classification factors.")

  if (lc$c1 == 0)
    stop("Some factor levels' combinations have zero frequency.")
  
  if (lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")
  
  if (lc$c3 == 0)
    stop("Some factor levels' combinations have additional replications.")

  if (lc$c4 == 1)
    stop("The data set is balanced. There are no missing values to estimate.")

  if (lc$pmis > maxp)
    stop(paste0("Too many missing values (", format(lc$pmis * 100, digits = 3), "%)."))
  
  if (sum(lc$nl < 2) > 0)
    stop("There are factors with only one level. This is is not a factorial experiment.")
  
  # Number of factors
  
  nf <- length(factors)

  # Get a copy of trait for estimated values and for temporary values
  
  trait.est <- paste0(trait, ".est")
  data[, trait.est] <- data[, trait]
  data[, "ytemp"] <- data[, trait]
  
  # Create expression for list of factors
  
  lf.expr <- 'list(data[, factors[1]]'
  
  for (i in 2:nf)
    lf.expr <- paste0(lf.expr, ', data[, factors[', i, ']]')
  
  lf.expr <- paste0(lf.expr, ')')
  
  # Compute means over replications

  tmeans <- tapply(data[, trait], eval(parse(text = lf.expr)), mean, na.rm = TRUE)

  # Store means in ytemp for missing values
  
  for (i in 1:length(data[, trait]))
    if (is.na(data[i, trait])) {
      expr <- paste0('tmeans[data[', i, ', factors[1]]')
      for (j in 2:nf)
        expr <- paste0(expr, ', data[', i, ', factors[', j, ']]')
      expr <- paste0(expr, ']')
      data[i, "ytemp"] <- eval(parse(text = expr))
    }
  
  # Estimate missing values for a crd
  
  if (design == "crd")
    data[, trait.est] <- data[, "ytemp"]
  
  # Estimate missing values for a rcbd
  
  if (design == "rcbd"){
    
    # Total number of treatments
    
    nt <- prod(lc$nl)
    
    # Vectors of estimated missing values to check convergence
    
    emv1 <- array(0, lc$nmis)
    emv2 <- array(0, lc$nmis)
    
    # Value to control convergence
    
    cc <- max(data[, trait], na.rm = TRUE)
    
    # Maximum number of iteration control
    
    cont <- 0
    
    # Iterate
    
    while (cc > max(data[, trait], na.rm = TRUE) * tol & cont < 100) {
      
      cont <- cont + 1
      
      for (i in 1:length(data[, trait]))
        if (is.na(data[i, trait])) {
      
          # Compute sums
          
          data[i, "ytemp"] <- NA
          sum1 <- tapply(data[, "ytemp"], eval(parse(text = lf.expr)), sum, na.rm = TRUE)
          sum2 <- tapply(data[, "ytemp"], data[, rep], sum, na.rm = TRUE)
          sum3 <- sum(data[, "ytemp"], na.rm = TRUE)
          
          # Get estimate
          
          expr <- paste0('sum1[data[', i, ', factors[1]]')
          for (j in 2:nf)
            expr <- paste0(expr, ', data[', i, ', factors[', j, ']]')
          expr <- paste0(expr, ']')
          
          mv.num <- nt * eval(parse(text = expr)) + lc$nr * sum2[data[i, rep]] - sum3
          mv.den <- nt * lc$nr - nt - lc$nr + 1

          data[i, trait.est] <- mv.num / mv.den
          
          data[i, "ytemp"] <- data[i, trait.est]
        }
      
      emv1 <- emv2
      emv2 <- data[is.na(data[, trait]), trait.est]
      cc <- max(abs(emv1 - emv2))
    }
  }
  
  # Return
  
  data[, c(factors, rep, trait, trait.est)]
  
}
