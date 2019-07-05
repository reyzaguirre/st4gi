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
#' @param nrs Number of rows to span the row of the plot.
#' @param method The method to fit the values. See details.
#' @param ind Logical. See details.
#' @param p The proportion of the check values differences used for the adjustment.
#' See details.
#' @param dfr The name of the data frame.
#' @details The values of the selected \code{trait} are adjusted using some mean
#' of the values of all the checks located on the row of the plot plus the \code{nrs}
#' rows at each side of the row of the plot. If \code{method = "flat"} the simple mean
#' of the checks is used for the adjustment, if \code{method = "weighted"} a weighted mean
#' is used for the adjustmen, where the checks closer to the plot get more weight.
#' If \code{p = 1} then the values are adjusted in the same proportion that the checks
#' vary around the field, for values lower than 1 the values are adjusted based on that
#' proportion over the checks variation, if \code{p = 0} then there is no adjustment.
#' If \code{ind = TRUE}, each check is centered around its own mean before the adjustment,
#' if \code{ind = FALSE}, then both checks are centered around the overall mean.
#' If \code{method = "flat"}, \code{ncb = 5}, \code{nrs = 1}, and \code{ind = FALSE}, then
#' it corresponds to the method proposed by Westcott.
#' If not specified, \code{nrs = floor(ncb / 2)}.
#' If the layout does not correspond with the Westcott method, then the observed values
#' are adjusted with the values of the checks planted nearby (this is in the rectangular
#' region definied by \code{floor(ncb / 2)} columns and rows at each side of
#' the plot) and a warning is issued.
#' @return It returns the adjusted values.
#' @author Raul Eyzaguirre.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @examples 
#' # Create design
#' dfr <- cr.w(1:1000, "Dag", "Cem", 40, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(dim(dfr)[1])
#' # Run adjustment
#' aj.w("y", "geno", "Dag", "Cem", "row", "col", dfr = dfr)
#' @export

aj.w <- function(trait, geno, ck1, ck2, row, col, ncb = 10, nrs = NULL,
                 method = c("weighted", "flat"), ind = TRUE, p = 0.5, dfr) {
  
  # Match arguments
  
  method <- match.arg(method)
  
  # Define nrs
  
  if (is.null(nrs))
    nrs <- floor(ncb / 2)
  
  # Error and warning messages
  
  out <- ck.pos(row, col, NULL, dfr)
  
  if (out$nplot > 0)
    stop("More than one genotype in the same position. Run check.pos to look over.")
  
  out <- ck.w(trait, geno, ck1, ck2, row, col, ncb, dfr)
  eval.cond <- sum(out$c1, out$c2, out$c3, out$c4, out$c5)
  
  if (out$c1 == 1)
    warning("There are plots in the columns of checks with other genotypes planted.")
  
  if (out$c2 == 1)
    warning("There are plots in the columns of genotypes with checks planted.")
  
  if (out$c3 == 1)
    warning("There are columns of checks with missing plots.")

  if (out$c4 == 1)
    warning("There are columns of checks without alternating checks.")
  
  if (out$c5 == 1)
    warning("There are plots with genotypes without a check plot to the left or to the right.")
  
  if (eval.cond > 0)
    warning("Adjusted values are obtained with the values of the checks nearby.")

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
  
  # Center check values and compute relative variation adjusted by p
  
  cond1 <- dfr[, geno] == ck1
  cond2 <- dfr[, geno] == ck2

  dfr[cond1, trait.aj] <- (dfr[cond1, trait.aj] - ck1.mean) / ck1.mean * p
  dfr[cond2, trait.aj] <- (dfr[cond2, trait.aj] - ck2.mean) / ck2.mean * p
  
  # Replace missing values with 0 (this is the centered mean)
  
  dfr[dfr[, geno] %in% c(ck1, ck2) & is.na(dfr[, trait.aj]), trait.aj] <- 0
  
  # Choose function for adjustment
  
  if (eval.cond == 0 & method == "weighted")
    foo <- foo.weig
  if (eval.cond == 0 & method == "flat")
    foo <- foo.flat
  if (eval.cond > 0) {
    nrs <- floor(ncb / 2)
    ncb <- floor(ncb / 2)
    foo <- foo.flat
  }
  
  # Make adjustment
  
  for (i in 1:dim(dfr)[1]) {
    if (dfr[i, geno] %in% c(ck1, ck2)) {
      af <- 0
    } else {
        af <- foo(i, trait.aj, geno, ck1, ck2, row, col, ncb, nrs, dfr)
    }
    dfr[i, trait.aj] <- dfr[i, trait.aj] / (1 + af)
  }
  
  # Return
  
  dfr
  
}

# A function for Westcott adjustment with method 2

foo.weig <- function(x, trait.aj, geno, ck1, ck2, row, col, ncb, nrs, dfr) {
  
  # Identify row and column for plot
  
  geno.row <- dfr[x, row]
  geno.col <- dfr[x, col]
  
  # Identify columns for checks (left and right)
  
  columns <- (geno.col - ncb):(geno.col + ncb)
  temp <- dfr[dfr[, row] == geno.row & dfr[, col] %in% columns & dfr[, geno] %in% c(ck1, ck2), ]
  col.lf <- min(temp[, col])
  col.rg <- max(temp[, col])
  
  # Identify rows for checks and define weights
  
  row.ck <- (geno.row - nrs):(geno.row + nrs)
  row.wg <- c(1:nrs, nrs + 1, nrs:1)
  
  # Delete nonexistent rows
  
  temp <- dfr[dfr[, col] == geno.col, ]
  valid.rows <- row.ck %in% temp[, row]
  row.ck <- row.ck[valid.rows]
  row.wg <- row.wg[valid.rows]
  
  # Check values on the left and right
  
  ck.lf <- dfr[dfr[, row] %in% row.ck & dfr[, col] == col.lf, c(row, trait.aj)]
  ck.rg <- dfr[dfr[, row] %in% row.ck & dfr[, col] == col.rg, c(row, trait.aj)]
  
  # Sort by row
  
  ck.lf <- ck.lf[sort.int(ck.lf[, row], index.return = TRUE)$ix, ]
  ck.rg <- ck.rg[sort.int(ck.rg[, row], index.return = TRUE)$ix, ]
  
  # Get left and right means for adjustment
  
  m.ck.lf <- sum(ck.lf[, trait.aj] * row.wg) / sum(row.wg)
  m.ck.rg <- sum(ck.rg[, trait.aj] * row.wg) / sum(row.wg)
  
  # Get adjustment factor
  
  af <- m.ck.lf + (geno.col - col.lf) * (m.ck.rg - m.ck.lf) / (col.rg - col.lf)
  
  # Return
  
  af
  
}

# A function for Westcott adjustment with method 1 or
# when layout does not follow the Westcott method

foo.flat <- function(x, trait.aj, geno, ck1, ck2, row, col, ncb, nrs, dfr) {
  
  # Identify row and column for plot
  
  geno.row <- dfr[x, row]
  geno.col <- dfr[x, col]
  
  # Identify all check values in the neighbourhood
  
  rows <- (geno.row - nrs):(geno.row + nrs)
  columns <- (geno.col - ncb):(geno.col + ncb)

  temp <- dfr[dfr[, row] %in% rows & dfr[, col] %in% columns & dfr[, geno] %in% c(ck1, ck2), ]

  # Get adjustment factor
  
  af <- mean(temp[, trait.aj])
  
  # Return
  
  af
  
}
