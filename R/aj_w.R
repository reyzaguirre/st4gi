#' Adjust values following the method of Westcott
#'
#' This function adjust the observed values of an experiment planted following
#' the method described by Westcott (1981) with a grid of checks.
#' @param trait The trait to adjust.
#' @param geno The genotypes.
#' @param ck1 Name of check 1.
#' @param ck2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param nrs Number of rows spanning the row of the plot, 3 or 5.
#' @param method The method to fit the values. See details.
#' @param p The proportion of the check values differences used for the adjustment.
#' See details.
#' @param ind Logical. If TRUE each check is centered around its own mean.
#' If FALSE both checks are centered around the overall mean.
#' @param dfr The name of the data frame.
#' @details The values of the selected \code{trait} are adjusted following
#' different methods. For \code{method = 1} each plot is adjusted with the
#' mean of the six or ten checks located in the three (if \code{nrs = 3}) or
#' five (if \code{nrs = 5}) rows spanning it (with \code{nrs = 3} it corresponds
#' tp the method proposed by Westcott). For \code{method = 2} each plot is
#' adjusted with the value of his position on the surface fitted using the six
#' or then checks located in the three or five rows spanning it, so the closest
#' checks receive more weight.
#' 
#' If \code{p = 1} then the values are adjusted in the same
#' proportion that the checks vary around the field. For values lower than 1
#' the values are adjusted based on that proportion over the checks variation.
#' If \code{p = 0} then there is no adjustment.
#' 
#' If the layout does not correspond with the Westcott method, then the observed values
#' are adjusted with the values of the checks planted nearby and a warning is issued.
#' @return It returns the adjusted values.
#' @author Raul Eyzaguirre.
#' @examples 
#' # Create design
#' dfr <- cr.w(1:1000, "Dag", "Cem", 40, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(1134)
#' # Run adjustment
#' aj.w("y", "geno", "Dag", "Cem", "row", "col", dfr = dfr)
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @export

aj.w <- function(trait, geno, ck1, ck2, row, col, nrs = 5, ncb = 10, method = 2,
                 p = 0.5, ind = TRUE, dfr) {
  
  # Error and warning messages
  
  out <- ck.pos(row, col, NULL, dfr)
  
  if (out$nplot > 0)
    stop("More than one genotype in the same position. Run check.pos to look over.")
  
  out <- ck.w(trait, geno, ck1, ck2, row, col, ncb, dfr)
  
  if (out$c1 == 1)
    warning("There are plots in the columns of checks with other genotypes planted.")
  
  if (out$c2 == 1)
    warning("There are plots in the columns of genotypes with checks planted.")
  
  if (out$c3 == 1)
    warning("There are columns of checks without alternating checks.")
  
  if (out$c4 == 1)
    warning("There are plots with genotypes without a check plot to the left or to the right.")
  
  if (out$c1 == 1 | out$c2 == 1 | out$c3 == 1 | out$c4 == 1)
    warning("Adjusted values are obtained with the values of the checks nearby.")

  # Save column names
  
  col.names <- colnames(dfr)
    
  # Get a copy of trait for the adjusted values
  
  trait.aj <- paste(trait, 'aj', sep = '.') 
  dfr[, trait.aj] <- dfr[, trait]

  # Compute means for checks
  
  if (ind) {
    ck1.mean <- mean(dfr[dfr[, geno] == ck1, trait], na.rm = TRUE)
    ck2.mean <- mean(dfr[dfr[, geno] == ck2, trait], na.rm = TRUE)
  } else {
    ck1.mean <- mean(dfr[dfr[, geno] %in% c(ck1, ck2), trait], na.rm = TRUE)
    ck2.mean <- mean(dfr[dfr[, geno] %in% c(ck1, ck2), trait], na.rm = TRUE)
  }
  
  # Center check values and compute %
  
  cond1 <- dfr[, geno] == ck1
  cond2 <- dfr[, geno] == ck2

  dfr[cond1, trait.aj] <- (dfr[cond1, trait.aj] - ck1.mean) / ck1.mean * p
  dfr[cond2, trait.aj] <- (dfr[cond2, trait.aj] - ck2.mean) / ck2.mean * p
  
  # Replace missing values with 0 (this is the centered mean)
  
  dfr[dfr[, geno] %in% c(ck1, ck2) & is.na(dfr[, trait.aj]), trait.aj] <- 0
  
  # Run Westcott adjustment or modify adjustment
  
  if (out$c1 == 0 & out$c2 == 0 & out$c3 == 0 & out$c4 == 0) {
    
    # Create columns for check centered values
    
    dfr[, ck1] <- NA
    dfr[, ck2] <- NA
    
    # Create columns for prior and posterior check centered values
    
    ck1.pri.1 <- paste(ck1, 'pri.1', sep = '.')
    dfr[, ck1.pri.1] <- NA
    ck1.pri.2 <- paste(ck1, 'pri.2', sep = '.')
    dfr[, ck1.pri.2] <- NA
    ck1.pos.1 <- paste(ck1, 'pos.1', sep = '.')
    dfr[, ck1.pos.1] <- NA
    ck1.pos.2 <- paste(ck1, 'pos.2', sep = '.')
    dfr[, ck1.pos.2] <- NA
    ck2.pri.1 <- paste(ck2, 'pri.1', sep = '.')
    dfr[, ck2.pri.1] <- NA
    ck2.pri.2 <- paste(ck2, 'pri.2', sep = '.')
    dfr[, ck2.pri.2] <- NA
    ck2.pos.1 <- paste(ck2, 'pos.1', sep = '.')
    dfr[, ck2.pos.1] <- NA
    ck2.pos.2 <- paste(ck2, 'pos.2', sep = '.')
    dfr[, ck2.pos.2] <- NA
    
    # Create columns for weigths
    
    ck1.w <- paste(ck1, 'w', sep = '.')
    dfr[, ck1.w] <- NA
    ck2.w <- paste(ck2, 'w', sep = '.')
    dfr[, ck2.w] <- NA
    
    # Arrange check values and weights
    
    for (i in 1:dim(dfr)[1]) {
      
      geno.row <- dfr[i, row]
      geno.col <- dfr[i, col]
      columns <- (geno.col - ncb):(geno.col + ncb)
      
      cond1 <- dfr[, col] %in% columns
      cond2 <- dfr[, geno] %in% c(ck1, ck2)
      
      temp <- dfr[dfr[, row] == geno.row & cond1 & cond2, c(geno, trait.aj, col)]
      
      if (dim(temp)[1] == 2) {
        
        # Checks on the row
        
        dfr[i, ck1] <- temp[temp[, geno] == ck1, trait.aj]
        dfr[i, ck2] <- temp[temp[, geno] == ck2, trait.aj]
        
        # Checks on row -2
        
        temp.pri <- dfr[dfr[, row] == geno.row - 2 & cond1 & cond2, c(geno, trait.aj, col)]
        
        if (dim(temp.pri)[1] == 2) {
          dfr[i, ck1.pri.2] <- temp.pri[temp.pri[, geno] == ck1, trait.aj]
          dfr[i, ck2.pri.2] <- temp.pri[temp.pri[, geno] == ck2, trait.aj]
        }
        
        # Checks on row -1
        
        temp.pri <- dfr[dfr[, row] == geno.row - 1 & cond1 & cond2, c(geno, trait.aj, col)]
        
        if (dim(temp.pri)[1] == 2) {
          dfr[i, ck1.pri.1] <- temp.pri[temp.pri[, geno] == ck2, trait.aj]
          dfr[i, ck2.pri.1] <- temp.pri[temp.pri[, geno] == ck1, trait.aj]
        }
        
        # Checks on row +1
        
        temp.pos <- dfr[dfr[, row] == geno.row + 1 & cond1 & cond2, c(geno, trait.aj, col)]
        
        if (dim(temp.pos)[1] == 2) {
          dfr[i, ck1.pos.1] <- temp.pos[temp.pos[, geno] == ck2, trait.aj]
          dfr[i, ck2.pos.1] <- temp.pos[temp.pos[, geno] == ck1, trait.aj]
        }
        
        # Checks on row +2  
        
        temp.pos <- dfr[dfr[, row] == geno.row + 2 & cond1 & cond2, c(geno, trait.aj, col)]
        
        if (dim(temp.pos)[1] == 2) {
          dfr[i, ck1.pos.2] <- temp.pos[temp.pos[, geno] == ck1, trait.aj]
          dfr[i, ck2.pos.2] <- temp.pos[temp.pos[, geno] == ck2, trait.aj]
        }
        
        # Weights for closest checks
        
        dfr[i, ck1.w] <- ncb + 1 - abs(temp[temp[, geno] == ck1, col] - geno.col)
        dfr[i, ck2.w] <- ncb + 1 - abs(temp[temp[, geno] == ck2, col] - geno.col)        
      }
    }
    
    # Adjust values with method 1
    
    if (method == 1) {
      chs <- c(ck1, ck2, ck1.pri.1, ck2.pri.1, ck1.pos.1, ck2.pos.1)
      if (nrs == 5)
        chs <- c(chs, ck1.pri.2, ck2.pri.2, ck1.pos.2, ck2.pos.2)
      af <- apply(dfr[, chs], 1, mean, na.rm = TRUE)
    }
    
    # Adjust values with method 2
    
    if (method == 2) {
      
      if (nrs == 3) {
        ck1.m.pri <- apply(dfr[, c(ck1, ck1.pri.1)], 1, mean, na.rm = TRUE)
        ck2.m.pri <- apply(dfr[, c(ck2, ck2.pri.1)], 1, mean, na.rm = TRUE)
        ck1.m.pos <- apply(dfr[, c(ck1, ck1.pos.1)], 1, mean, na.rm = TRUE)
        ck2.m.pos <- apply(dfr[, c(ck2, ck2.pos.1)], 1, mean, na.rm = TRUE)
      }
      
      if (nrs == 5) {
        foo <- function(x) x * c(1.5, 2, 1)
        ck1.w.pri <- t(apply(!is.na(dfr[, c(ck1, ck1.pri.1, ck1.pri.2)]), 1, foo))
        ck2.w.pri <- t(apply(!is.na(dfr[, c(ck2, ck2.pri.1, ck2.pri.2)]), 1, foo))
        ck1.w.pos <- t(apply(!is.na(dfr[, c(ck1, ck1.pos.1, ck1.pos.2)]), 1, foo))
        ck2.w.pos <- t(apply(!is.na(dfr[, c(ck2, ck2.pos.1, ck2.pos.2)]), 1, foo))
        
        ck1.m.pri <- dfr[, c(ck1, ck1.pri.1, ck1.pri.2)] * ck1.w.pri
        ck2.m.pri <- dfr[, c(ck2, ck2.pri.1, ck2.pri.2)] * ck2.w.pri
        ck1.m.pos <- dfr[, c(ck1, ck1.pos.1, ck1.pos.2)] * ck1.w.pos
        ck2.m.pos <- dfr[, c(ck2, ck2.pos.1, ck2.pos.2)] * ck2.w.pos
        
        ck1.m.pri <- apply(ck1.m.pri, 1, sum, na.rm = TRUE)
        ck2.m.pri <- apply(ck2.m.pri, 1, sum, na.rm = TRUE)
        ck1.m.pos <- apply(ck1.m.pos, 1, sum, na.rm = TRUE)
        ck2.m.pos <- apply(ck2.m.pos, 1, sum, na.rm = TRUE)
        
        ck1.w.pri <- apply(ck1.w.pri, 1, sum, na.rm = TRUE)
        ck2.w.pri <- apply(ck2.w.pri, 1, sum, na.rm = TRUE)
        ck1.w.pos <- apply(ck1.w.pos, 1, sum, na.rm = TRUE)
        ck2.w.pos <- apply(ck2.w.pos, 1, sum, na.rm = TRUE)
        
        ck1.m.pri <- ck1.m.pri / ck1.w.pri
        ck2.m.pri <- ck2.m.pri / ck2.w.pri
        ck1.m.pos <- ck1.m.pos / ck1.w.pos
        ck2.m.pos <- ck2.m.pos / ck2.w.pos
      }
      
      l.pri <- (ck1.m.pri * dfr[, ck1.w] + ck2.m.pri * dfr[, ck2.w]) / (dfr[, ck1.w] + dfr[, ck2.w])
      l.pos <- (ck1.m.pos * dfr[, ck1.w] + ck2.m.pos * dfr[, ck2.w]) / (dfr[, ck1.w] + dfr[, ck2.w])
      
      af <- apply(cbind(l.pri, l.pos), 1, mean)
    }
    
    # Make adjustment
    
    dfr[, trait.aj] <- dfr[, trait.aj] / (1 + af)
    
  } else { 

    for (i in 1:dim(dfr)[1]) {
      
      geno.row <- dfr[i, row]
      geno.col <- dfr[i, col]
      rows <- (geno.row - nrs %/% 2):(geno.row + nrs %/% 2)
      columns <- (geno.col - ncb):(geno.col + ncb)
      
      cond1 <- dfr[, col] %in% columns & dfr[, row] %in% rows
      cond2 <- dfr[, geno] %in% c(ck1, ck2)
      
      temp <- dfr[cond1 & cond2, trait.aj]
      
      if (length(temp) > 0 & !(dfr[i, geno] %in% c(ck1, ck2))) {
        af <- mean(temp)
        dfr[i, trait.aj] <- dfr[i, trait.aj] / (1 + af)
      }
    }
  }
  
  # Return
  
  dfr[, unique(c(col.names, trait.aj))]
  
}
