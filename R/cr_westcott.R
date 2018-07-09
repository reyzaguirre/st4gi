#' Unreplicated experiment with a grid of checks
#'
#' This function creates the fieldbook and fieldplan for an unreplicated
#' experiment with genotypes randomly allocated on a field with checks
#' following the method described by Westcott (1981).
#' @param geno The list of genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param nc Number of columns.
#' @param ncb Number of columns between two check columns (default is 10).
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field between equally spaced
#' columns of two alternating check varieties. Check columns are planted each
#' \code{ncb} columns. The specified total number of columns \code{nc} is the maximum
#' available number on the field, the actual number could be less. 
#' @return It returns the fieldbook and fieldplan.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @examples
#' cr.w(1:100, "A", "B", 100)
#' @export

cr.w <- function(geno, ch1, ch2, nc, ncb = 10) {
  
  # As character
  
  geno <- as.character(geno)

  # Error messages
  
  ng <- length(geno)

  if (ng < ncb)
    stop(paste("Include at least", ncb, "genotypes."))
  
  # Dimensions
  
  nb <- ceiling(ng / ncb)            # Number of blocks
  nbr <- floor((nc - 1) / (ncb + 1)) # Number of blocks per row
  nc <- nbr * (ncb + 1) + 1          # Actual number of columns
  nr <- ceiling(ng / nbr / ncb)      # Number of rows
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Include checks, selected columns
  
  plan[seq(1, nr, 2), seq(1, nc, 2 + 2 * ncb)] <- ch1
  plan[seq(1, nr, 2), seq(2 + ncb, nc, 2 + 2 * ncb)] <- ch2
  
  if (nr > 1) {
    plan[seq(2, nr, 2), seq(1, nc, 2 + 2 * ncb)] <- ch2
    plan[seq(2, nr, 2), seq(2 + ncb, nc, 2 + 2 * ncb)] <- ch1
  }
  
  # Include genotypes at random
  
  geno <- sample(geno)
  
  k <- 1
  
  for (i in 1:nr)
    for (b in 1:nbr)
      for (j in (b + 1 + (b - 1) * ncb):(b * (1 + ncb))) {
        plan[i, j] <- geno[k]
        k <- k + 1
      }

  # Delete non necessary checks
  
  neb <- nbr * nr - nb  # Number of empty blocks
  if (neb > 0)
    plan[nr, (nc - (ncb + 1) * neb + 1):nc] <- NA

  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  book <- data.frame(plot = NA, row, col, geno = c(t(plan)), stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]
  book$plot <- 1:dim(book)[1]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
}
