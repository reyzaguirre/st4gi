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
#' @param nr Number of rows used to fit values, 3 or 5.
#' @param method The method to fit the values. See details.
#' @param p The proportion of the check values differences used for the adjustment.
#' See details.
#' @param ind Logical. If TRUE each check is centered around its own mean.
#' If FALSE both checks are centered around the overall mean.
#' @param data The name of the data frame.
#' @author Raul Eyzaguirre.
#' @details The values of the selected \code{trait} are adjusted following
#' different methods. For \code{method = 1} each plot is adjusted with the
#' mean of the six or ten checks located in the three (if \code{nr = 3}) or
#' five (if \code{nr = 5}) rows spanning it (with \code{nr = 3} it corresponds
#' tp the method proposed by Westcott). For \code{method = 2} each plot is
#' adjusted with the value of his position on the surface fitted using the six
#' or then checks located in the three or five rows spanning it, so the closest
#' checks receive more weight.
#' 
#' If \code{p = 1} then the values are adjusted in the same
#' proportion that the checks vary around the field. For values lower than 1
#' the values are adjusted based on that proportion over the checks variation.
#' If \code{p = 0} then there is no adjustment.
#' @return It returns the adjusted values.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @export

aj.wd <- function(trait, geno, ch1, ch2, row, col, nr = 5, ncb = 10, method = 2,
                  p = 0.5, ind = TRUE, data) {
  
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
  
  cond1 <- data[, geno] == ch1
  cond2 <- data[, geno] == ch2

  data[cond1, trait.aj] <- (data[cond1, trait.aj] - ch1.mean) / ch1.mean * p
  data[cond2, trait.aj] <- (data[cond2, trait.aj] - ch2.mean) / ch2.mean * p
  
  # Replace missing values with 0 (this is the centered mean)
  
  data[data[, geno] %in% c(ch1, ch2) & is.na(data[, trait.aj]), trait.aj] <- 0
  
  # Create columns for check centered values
  
  data[, ch1] <- NA
  data[, ch2] <- NA

  # Create columns for prior and posterior check centered values
    
  ch1.pri.1 <- paste(ch1, 'pri.1', sep = '.')
  data[, ch1.pri.1] <- NA
  ch1.pri.2 <- paste(ch1, 'pri.2', sep = '.')
  data[, ch1.pri.2] <- NA
  ch1.pos.1 <- paste(ch1, 'pos.1', sep = '.')
  data[, ch1.pos.1] <- NA
  ch1.pos.2 <- paste(ch1, 'pos.2', sep = '.')
  data[, ch1.pos.2] <- NA
  ch2.pri.1 <- paste(ch2, 'pri.1', sep = '.')
  data[, ch2.pri.1] <- NA
  ch2.pri.2 <- paste(ch2, 'pri.2', sep = '.')
  data[, ch2.pri.2] <- NA
  ch2.pos.1 <- paste(ch2, 'pos.1', sep = '.')
  data[, ch2.pos.1] <- NA
  ch2.pos.2 <- paste(ch2, 'pos.2', sep = '.')
  data[, ch2.pos.2] <- NA
  
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
    
    cond1 <- data[, col] %in% columns
    cond2 <- data[, geno] %in% c(ch1, ch2)
      
    temp <- data[data[, row] == geno.row & cond1 & cond2, c(geno, trait.aj, col)]
      
    if (dim(temp)[1] == 2) {
      
      # Checks on the row
      
      data[i, ch1] <- temp[temp[, geno] == ch1, trait.aj]
      data[i, ch2] <- temp[temp[, geno] == ch2, trait.aj]
      
      # Checks on row -2
      
      temp.pri <- data[data[, row] == geno.row - 2 & cond1 & cond2, c(geno, trait.aj, col)]
        
      if (dim(temp.pri)[1] == 2) {
        data[i, ch1.pri.2] <- temp.pri[temp.pri[, geno] == ch1, trait.aj]
        data[i, ch2.pri.2] <- temp.pri[temp.pri[, geno] == ch2, trait.aj]
      }
      
      # Checks on row -1
      
      temp.pri <- data[data[, row] == geno.row - 1 & cond1 & cond2, c(geno, trait.aj, col)]
      
      if (dim(temp.pri)[1] == 2) {
        data[i, ch1.pri.1] <- temp.pri[temp.pri[, geno] == ch2, trait.aj]
        data[i, ch2.pri.1] <- temp.pri[temp.pri[, geno] == ch1, trait.aj]
      }
      
      # Checks on row +1

      temp.pos <- data[data[, row] == geno.row + 1 & cond1 & cond2, c(geno, trait.aj, col)]
        
      if (dim(temp.pos)[1] == 2) {
        data[i, ch1.pos.1] <- temp.pos[temp.pos[, geno] == ch2, trait.aj]
        data[i, ch2.pos.1] <- temp.pos[temp.pos[, geno] == ch1, trait.aj]
      }
      
      # Checks on row +2  

      temp.pos <- data[data[, row] == geno.row + 2 & cond1 & cond2, c(geno, trait.aj, col)]
      
      if (dim(temp.pos)[1] == 2) {
        data[i, ch1.pos.2] <- temp.pos[temp.pos[, geno] == ch1, trait.aj]
        data[i, ch2.pos.2] <- temp.pos[temp.pos[, geno] == ch2, trait.aj]
      }
      
      # Weights for closest checks

      data[i, ch1.w] <- ncb + 1 - abs(temp[temp[, geno] == ch1, col] - geno.col)
      data[i, ch2.w] <- ncb + 1 - abs(temp[temp[, geno] == ch2, col] - geno.col)        
    }
  }
    
  # Adjust values with method 1
    
  if (method == 1) {
    chs <- c(ch1, ch2, ch1.pri.1, ch2.pri.1, ch1.pos.1, ch2.pos.1)
    if (nr == 5)
      chs <- c(chs, ch1.pri.2, ch2.pri.2, ch1.pos.2, ch2.pos.2)
    af <- apply(data[, chs], 1, mean, na.rm = TRUE)
  }
    
  # Adjust values with method 2
      
  if (method == 2) {
    
    if (nr == 3) {
      ch1.m.pri <- apply(data[, c(ch1, ch1.pri.1)], 1, mean, na.rm = TRUE)
      ch2.m.pri <- apply(data[, c(ch2, ch2.pri.1)], 1, mean, na.rm = TRUE)
      ch1.m.pos <- apply(data[, c(ch1, ch1.pos.1)], 1, mean, na.rm = TRUE)
      ch2.m.pos <- apply(data[, c(ch2, ch2.pos.1)], 1, mean, na.rm = TRUE)
    }
      
    if (nr == 5) {
      foo <- function(x) x * c(1.5, 2, 1)
      ch1.w.pri <- t(apply(!is.na(data[, c(ch1, ch1.pri.1, ch1.pri.2)]), 1, foo))
      ch2.w.pri <- t(apply(!is.na(data[, c(ch2, ch2.pri.1, ch2.pri.2)]), 1, foo))
      ch1.w.pos <- t(apply(!is.na(data[, c(ch1, ch1.pos.1, ch1.pos.2)]), 1, foo))
      ch2.w.pos <- t(apply(!is.na(data[, c(ch2, ch2.pos.1, ch2.pos.2)]), 1, foo))
      
      ch1.m.pri <- data[, c(ch1, ch1.pri.1, ch1.pri.2)] * ch1.w.pri
      ch2.m.pri <- data[, c(ch2, ch2.pri.1, ch2.pri.2)] * ch2.w.pri
      ch1.m.pos <- data[, c(ch1, ch1.pos.1, ch1.pos.2)] * ch1.w.pos
      ch2.m.pos <- data[, c(ch2, ch2.pos.1, ch2.pos.2)] * ch2.w.pos
      
      ch1.m.pri <- apply(ch1.m.pri, 1, sum, na.rm = TRUE)
      ch2.m.pri <- apply(ch2.m.pri, 1, sum, na.rm = TRUE)
      ch1.m.pos <- apply(ch1.m.pos, 1, sum, na.rm = TRUE)
      ch2.m.pos <- apply(ch2.m.pos, 1, sum, na.rm = TRUE)
      
      ch1.w.pri <- apply(ch1.w.pri, 1, sum, na.rm = TRUE)
      ch2.w.pri <- apply(ch2.w.pri, 1, sum, na.rm = TRUE)
      ch1.w.pos <- apply(ch1.w.pos, 1, sum, na.rm = TRUE)
      ch2.w.pos <- apply(ch2.w.pos, 1, sum, na.rm = TRUE)
      
      ch1.m.pri <- ch1.m.pri / ch1.w.pri
      ch2.m.pri <- ch2.m.pri / ch2.w.pri
      ch1.m.pos <- ch1.m.pos / ch1.w.pos
      ch2.m.pos <- ch2.m.pos / ch2.w.pos
    }

    l.pri <- (ch1.m.pri * data[, ch1.w] + ch2.m.pri * data[, ch2.w]) / (data[, ch1.w] + data[, ch2.w])
    l.pos <- (ch1.m.pos * data[, ch1.w] + ch2.m.pos * data[, ch2.w]) / (data[, ch1.w] + data[, ch2.w])
      
    af <- apply(cbind(l.pri, l.pos), 1, mean)
  }

  # Make adjustment
  
  data[, trait.aj] <- data[, trait.aj] / (1 + af)
  
  # Return
  
  data[, unique(c(col.names, trait.aj))]
  
}
