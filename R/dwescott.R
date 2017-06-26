#' Unreplicated design with a grid of checks
#'
#' This function creates the fieldbook and fieldplan for an unreplicated design
#' with genotypes randomly allocated on a field with checks following
#' the method described on Westcott (1981).
#' @param geno The list of genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param nr Number of rows.
#' @param ncb Number of columns between two check columns.
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

dw <- function(geno, ch1, ch2, nr, ncb) {
  
  # Check ncb is even
  
  if (ncb %% 2 == 1)
    stop("The number of columns must be even.")
  
  # Dimensions
  
  ng <- length(geno)           # Number of genotypes
  nb <- ceiling(ng / nr / ncb) # Number of blocks
  nc <- (ncb + 1) * nb + 1     # Number of columns
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))
  
  # Include checks, selected columns
  
  plan[seq(1, nr, 2), seq(1, nc, 2 + 2 * ncb)] <- ch1
  plan[seq(2, nr, 2), seq(1, nc, 2 + 2 * ncb)] <- ch2
  plan[seq(1, nr, 2), seq(2 + ncb, nc, 2 + 2 * ncb)] <- ch2
  plan[seq(2, nr, 2), seq(2 + ncb, nc, 2 + 2 * ncb)] <- ch1
  
  # Include genotypes at random
  
  geno <- sample(geno)
  
  k <- 1
  
  for (b in 1:nb)
    for (i in 1:nr)
      for (j in (b + 1 + (b - 1) * ncb):(b * (1 + ncb))) {
        if (is.na(plan[i, j])) {
          plan[i, j] <- geno[k]
          k <- k + 1
        }
      }
  
  # Delete non necessary checks
  
  nrc <- sum(!is.na(plan[, nc - 1])) + 1 # number of rows with clones in the last column
  if (nrc <= nr)
    plan[nrc:nr, nc] <- NA
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  book <- data.frame(plot = 1:(nr * nc),
                     row, col, geno = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
}
