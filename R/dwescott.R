#' Unreplicated design with a grid of checks
#'
#' This function creates the fieldbook and fieldplan for an unreplicated design
#' with genotypes randomly allocated on a field with checks following
#' the method described on Westcott (1981).
#' @param geno The list of genotypes.
#' @param ch1 Name of check 1.
#' @param ch2 Name of check 2.
#' @param nc Number of columns.
#' @param ncb Number of columns between two check columns (default is 10).
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field between equally spaced
#' columns of two alternating check varieties. Check columns are planted each
#' \code{ncb} columns where \code{ncb} must be even. The specified total number of
#' columns \code{nc} is the maximum available number, the actual number could be less. 
#' @return It returns the fieldbook and fieldplan.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @examples
#' dw(1:100, "A", "B", 100)
#' @export

dw <- function(geno, ch1, ch2, nc, ncb = 10) {
  
  # Check ncb is even
  
  if (ncb %% 2 == 1)
    stop("The number of columns must be even.")
  
  # Dimensions
  
  ng <- length(geno)                 # Number of genotypes
  nb <- ceiling(ng / ncb)            # Number of blocks
  nbr <- floor((nc - 1) / (ncb + 1)) # Number of blocks per row
  nc <- nbr * (ncb + 1) + 1          # Actual number of columns
  nr <- ceiling(ng / nbr / ncb)      # Number of rows
  
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

  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  book <- data.frame(plot = NA, row, col, geno = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  book$plot <- 1:dim(book)[1]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
}
