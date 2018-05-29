#' Adjust values with a grid of checks
#'
#' This function adjust the observed values of an experiment planted following
#' the method described on Westcott (1981) with a grid of checks.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param method The method to fit the values. See details.
#' @param w The weight from 0 to 1 given to the check values for the adjustment. See details.
#' @param ind Logical. If TRUE each check is centered around its own mean. If FALSE both
#' checks are centered around the overall mean.
#' @param data The name of the data frame.
#' @author Raul Eyzaguirre.
#' @details The values of the selected \code{trait} are adjusted following different
#' methods. For \code{method = 1} each plot is adjusted with the mean of the six checks
#' located in the three rows spanning it (this is the method proposed by Westcott).
#' For \code{method = 2} each plot is adjusted with the value of his position on the
#' surface fitted using the six checks located in the three rows spanning it, so the closest
#' checks receive more weight.
#' 
#' \code{w} gives the weight given to the checks for the adjustmen. If \code{w = 1} then the
#' values are adjusted in the same proportion that the checks vary around the field. For
#' values lower than 1 the values are adjusted based on that proportion over the checks variation.
#' @return It returns the adjusted values.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @export

aj.wd <- function(trait, geno, ch1, ch2, row, col, ncb, method = 1, w = 0.25,
                  ind = TRUE, data) {
  
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
  
  if (ind) {
    ch1.mean <- mean(data[data[, geno] == ch1, trait], na.rm = TRUE)
    ch2.mean <- mean(data[data[, geno] == ch2, trait], na.rm = TRUE)
  } else {
    ch1.mean <- mean(data[data[, geno] %in% c(ch1, ch2), trait], na.rm = TRUE)
    ch2.mean <- mean(data[data[, geno] %in% c(ch1, ch2), trait], na.rm = TRUE)
  }
  
  # Center check values and compute %
  
  data[data[, geno] == ch1, trait.aj] <- (data[data[, geno] == ch1, trait.aj] - ch1.mean) / ch1.mean * w
  data[data[, geno] == ch2, trait.aj] <- (data[data[, geno] == ch2, trait.aj] - ch2.mean) / ch2.mean * w
  
  # Replace missing values with 0 (this is the centered mean)
  
  data[data[, geno] %in% c(ch1, ch2) & is.na(data[, trait.aj]), trait.aj] <- 0
  
  # Create columns for check centered values
  
  data[, ch1] <- NA
  data[, ch2] <- NA

  # Create columns for prior and posterior check centered values
    
  ch1.pri <- paste(ch1, 'pri', sep = '.')
  data[, ch1.pri] <- NA
  ch1.pos <- paste(ch1, 'pos', sep = '.')
  data[, ch1.pos] <- NA
  ch2.pri <- paste(ch2, 'pri', sep = '.')
  data[, ch2.pri] <- NA
  ch2.pos <- paste(ch2, 'pos', sep = '.')
  data[, ch2.pos] <- NA
    
  # Create columns for weigths
    
  ch1.w <- paste(ch1, 'w', sep = '.')
  data[, ch1.w] <- NA
  ch2.w <- paste(ch2, 'w', sep = '.')
  data[, ch2.w] <- NA
    
  # Arrange check values and weights
    
  for(i in 1:dim(data)[1]) {
      
    geno.row <- data[i, row]
    geno.col <- data[i, col]
    columns <- (geno.col - ncb):(geno.col + ncb)
      
    temp <- data[data[, row] == geno.row & data[, col] %in% columns & data[, geno] %in% c(ch1, ch2), c(geno, trait.aj, col)]
      
    if (dim(temp)[1] == 2) {
        
      data[i, ch1] <- temp[temp[, geno] == ch1, trait.aj]
      data[i, ch2] <- temp[temp[, geno] == ch2, trait.aj]
        
      temp.pri <- data[data[, row] == geno.row - 1 & data[, col] %in% columns & data[, geno] %in% c(ch1, ch2), c(geno, trait.aj, col)]
        
      if (dim(temp.pri)[1] == 2) {
        data[i, ch1.pri] <- temp.pri[temp.pri[, geno] == ch2, trait.aj]
        data[i, ch2.pri] <- temp.pri[temp.pri[, geno] == ch1, trait.aj]
      }
        
      temp.pos <- data[data[, row] == geno.row + 1 & data[, col] %in% columns & data[, geno] %in% c(ch1, ch2), c(geno, trait.aj, col)]
        
      if (dim(temp.pos)[1] == 2) {
        data[i, ch1.pos] <- temp.pos[temp.pos[, geno] == ch2, trait.aj]
        data[i, ch2.pos] <- temp.pos[temp.pos[, geno] == ch1, trait.aj]
      }
        
      data[i, ch1.w] <- ncb + 1 - abs(temp[temp[, geno] == ch1, col] - geno.col)
      data[i, ch2.w] <- ncb + 1 - abs(temp[temp[, geno] == ch2, col] - geno.col)        
    }
  }
    
  # Adjust values with method 3
    
  if (method == 1)
    af <- apply(data[, c(ch1, ch2, ch1.pri, ch2.pri, ch1.pos, ch2.pos)], 1, mean, na.rm = TRUE)

  # Adjust values with method 4
      
  if (method == 2) {
    m.pri.1 <- apply(data[, c(ch1, ch1.pri)], 1, mean, na.rm = TRUE)
    m.pri.2 <- apply(data[, c(ch2, ch2.pri)], 1, mean, na.rm = TRUE)
    m.pos.1 <- apply(data[, c(ch1, ch1.pos)], 1, mean, na.rm = TRUE)
    m.pos.2 <- apply(data[, c(ch2, ch2.pos)], 1, mean, na.rm = TRUE)
      
    l.pri <- (m.pri.1 * data[, ch1.w] + m.pri.2 * data[, ch2.w]) / (data[, ch1.w] + data[, ch2.w])
    l.pos <- (m.pos.1 * data[, ch1.w] + m.pos.2 * data[, ch2.w]) / (data[, ch1.w] + data[, ch2.w])
      
    af <- apply(cbind(l.pri, l.pos), 1, mean)
  }

  # Make adjustment
  
  data[, trait.aj] <- data[, trait.aj] / (1 + af)
  
  # Return
  
  data[, unique(c(col.names, trait.aj))]
  
}
