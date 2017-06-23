#' Unreplicated design with a grid of checks
#'
#' This function creates the fieldbook and fieldplan for an unreplicated design
#' with genotypes randomly allocated on a field with checks following
#' the method described on Westcott (1981).
#' @param geno The list of genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param nr Number of rows.
#' @param nc Number of columns between two check columns.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field with equally spaced
#' columns of two alternating check varieties.
#' @return It returns the fieldbook and fieldplan.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @examples
#' dw(1:100, "A", "B", 5, 10)
#' @export

dw <- function(geno, ch1, ch2, nr, nc) {
  
  # Check nc is even
  
  if (nc %% 2 == 1)
    stop("The number of columns must be even.")
  
  # Dimentions
  
  ng <- length(geno)          # Number of genotypes
  nb <- ceiling(ng / nr / nc) # Number of blocks
  tnc <- (nc + 1) * nb + 1    # Total number of columns
  
  # fieldplan array
  
  fp <- array(dim = c(nr, tnc))
  
  # Include checks, selected columns
  
  fp[seq(1, nr, 2), seq(1, tnc, 2 + 2 * nc)] <- ch1
  fp[seq(2, nr, 2), seq(1, tnc, 2 + 2 * nc)] <- ch2
  fp[seq(1, nr, 2), seq(2 + nc, tnc, 2 + 2 * nc)] <- ch2
  fp[seq(2, nr, 2), seq(2 + nc, tnc, 2 + 2 * nc)] <- ch1
  
  # Include genotypes at random
  
  geno <- sample(geno)
  k <- 1
  
  for (b in 1:nb)
    for (i in 1:nr)
      for (j in (b + 1 + (b - 1) * nc):(b * (1 + nc))) {
        if (is.na(fp[i, j])) {
          fp[i, j] <- geno[k]
          k <- k + 1
        }
      }
  
  # Delete non necessary checks
  
  nrc <- sum(!is.na(fp[, tnc - 1])) + 1 # number of rows with clones in the last column
  if (nrc <= nr)
    fp[nrc:nr, tnc] <- NA
  
  # row and column names
  
  rownames(fp) <- paste("row", 1:nr)
  colnames(fp) <- paste("col", 1:tnc)

  # Create fielbook
  
  row <- as.integer(gl(dim(fp)[1], dim(fp)[2]))
  col <- rep(1:dim(fp)[2], dim(fp)[1])
  book <- data.frame(plot = 1:(dim(fp)[1] * dim(fp)[2]),
                     row, col, geno = c(t(fp)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(fp = fp, book = book)
}
