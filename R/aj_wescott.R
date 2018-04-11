#' Adjust values with a grid of checks
#'
#' This function adjust the observed values of an experiment planted following
#' the method described on Wescott (1981) with a grid of checks.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param method The method to fit the values. See details.
#' @param data The name of the data frame.
#' @author Raul Eyzaguirre.
#' @details The values of the selected \code{trait} are adjusted following different
#' methods. For \code{method = 1} each value is adjusted with the mean of the checks
#' located to the left and right.
#' @return It returns the adjusted values.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @export

aj.wd <- function(trait, geno, ch1, ch2, row, col, ncb, method = 1, data) {
  
  # Error messages
  
  out <- check.rc(row, col, data = data)
  
  if (out$nplot > 0)
    stop("More than one genotype in the same position. Run check.pos to see.")
  
  out <- check.wd(trait, geno, ch1, ch2, row, col, ncb, data)
  
  if (out$c1 == 1)
    stop("There are plots in the columns of checks with other genotypes planted.")
  
  if (out$c2 == 1)
    stop("The last column in the field does not have checks.")
  
  if (out$c3 == 1)
    stop("There are plots in the columns of genotypes with checks planted.")
  
  if (out$c4 == 1)
    stop("There are columns of checks without alternating checks.")
  
  if (out$c5 == 1)
    stop("There are plots with genotypes without a check plot to the left or to the right.")
  
  # Save column names
  
  col.names <- colnames(data)
    
  # Get a copy of trait for the adjusted values
  
  trait.aj <- paste(trait, 'aj', sep = '.') 
  data[, trait.aj] <- data[, trait]

  # Compute means for checks
  
  ch1.mean <- mean(data[data[, geno] == ch1, trait], na.rm = TRUE)
  ch2.mean <- mean(data[data[, geno] == ch2, trait], na.rm = TRUE)
  
  # Center check values
  
  data[data[, geno] == ch1, trait.aj] <- data[data[, geno] == ch1, trait.aj] - ch1.mean
  data[data[, geno] == ch2, trait.aj] <- data[data[, geno] == ch2, trait.aj] - ch2.mean
  
  # Replace missing values with 0 (this is the centered mean)
  
  data[data[, geno] %in% c(ch1, ch2) & is.na(data[, trait.aj]), trait.aj] <- 0
  
  # Create columns for check centered values
  
  data[, ch1] <- NA
  data[, ch2] <- NA

  # Adjust values for method 1
  
  if (method == 1) {

    # Arrange check values properly
    
    for(i in 1:dim(data)[1]) {
      geno.row <- data[i, row]
      columns <- (data[i, col] - ncb):(data[i, col] + ncb)
      temp <- data[data[, row] == geno.row & data[, col] %in% columns & data[, geno] %in% c(ch1, ch2), c(geno, trait.aj)]
      if (dim(temp)[1] == 2) {
        data[i, ch1] <- temp[temp[, geno] == ch1, trait.aj]
        data[i, ch2] <- temp[temp[, geno] == ch2, trait.aj]
      }
    }
    
    # Adjust values with the mean of the centered check values
    
    data[, trait.aj] <- data[, trait.aj] - (data[, ch1] + data[, ch2]) / 2
    
  }
  
  # Adjust values for method 2
  
  if (method == 2) {
    
    # Create columns for weigths

    ch1.w <- paste(ch1, 'w', sep = '.')
    data[, ch1.w] <- NA
    ch2.w <- paste(ch2, 'w', sep = '.')
    data[, ch2.w] <- NA
    
    # Arrange check values properly
    
    for(i in 1:dim(data)[1]) {
      geno.row <- data[i, row]
      geno.col <- data[i, col]
      columns <- (geno.col - ncb):(geno.col + ncb)
      temp <- data[data[, row] == geno.row & data[, col] %in% columns & data[, geno] %in% c(ch1, ch2), c(geno, trait.aj, col)]
      if (dim(temp)[1] == 2) {
        data[i, ch1] <- temp[temp[, geno] == ch1, trait.aj]
        data[i, ch2] <- temp[temp[, geno] == ch2, trait.aj]
        data[i, ch1.w] <- ncb + 1 - abs(temp[temp[, geno] == ch1, col] - geno.col)
        data[i, ch2.w] <- ncb + 1 - abs(temp[temp[, geno] == ch2, col] - geno.col)        
      }
    }
    
    # Adjust values with the linear interpolation of the centered check values
    
    data[, trait.aj] <- data[, trait.aj] - (data[, ch1] * data[, ch1.w] + data[, ch2] * data[, ch2.w]) / (data[, ch1.w] + data[, ch2.w])
    
  }
  
  # Return
  
  data[, unique(c(col.names, trait.aj))]
  
}

