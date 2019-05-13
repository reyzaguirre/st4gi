#' Unreplicated experiment with a grid of checks
#'
#' This function creates the fieldbook and fieldplan for an unreplicated
#' experiment with genotypes randomly allocated on a field with checks
#' following the method described by Westcott (1981).
#' @param geno The list of genotypes.
#' @param ck1 Name of check 1.
#' @param ck2 Name of check 2.
#' @param nc Number of available columns on the field.
#' @param ncb Number of columns between two check columns (default is 10).
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @details The genotypes are randomly allocated on a field between equally spaced
#' columns of two alternating check varieties. Check columns are planted each
#' \code{ncb} columns.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @references
#' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
#' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
#' INRA Poitier, France, pp 91-95.
#' @examples
#' cr.w(1:100, "A", "B", 23)
#' @export

cr.w <- function(geno, ck1, ck2, nc = NULL, ncb = 10, serpentine = c("yes", "no")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)

  # Error messages
  
  ng <- length(geno)

  if (ng < ncb)
    stop(paste("Include at least", ncb, "genotypes."))
  
  # Number of columns
  
  if (is.null(nc))
    nc <- max(ncb + 2, round(sqrt(ng)))
  
  # Dimensions
  
  nb <- ceiling(ng / ncb)            # Number of incomplete blocks
  nbr <- floor((nc - 1) / (ncb + 1)) # Number of incomplete blocks per row
  nc <- nbr * (ncb + 1) + 1          # Actual number of columns
  nr <- ceiling(ng / nbr / ncb)      # Number of rows
  
  # Fieldplan array
  
  plan.id <- fp(nr, nc, serpentine)
  
  plan <- array(dim = c(nr, nc))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
  # Include checks, selected columns
  
  plan[seq(1, nr, 2), seq(1, nc, 2 + 2 * ncb)] <- ck1
  plan[seq(1, nr, 2), seq(2 + ncb, nc, 2 + 2 * ncb)] <- ck2
  
  if (nr > 1) {
    plan[seq(2, nr, 2), seq(1, nc, 2 + 2 * ncb)] <- ck2
    plan[seq(2, nr, 2), seq(2 + ncb, nc, 2 + 2 * ncb)] <- ck1
  }
  
  # Random order for genotypes
  
  geno <- sample(geno)
  
  # Get positions for genotypes
  
  temp <- plan.id
  temp[!is.na(plan)] <- NA
  temp <- c(temp)
  temp <- sort(temp[!is.na(temp)])
  temp <- temp[1:ng]
  
  # Create a data frame with positions and genotypes
  
  temp <- data.frame(temp, geno, stringsAsFactors = FALSE)

  # Asign genotypes to field plan
    
  for (k in 1:ng) {
    plan[plan.id == temp[k, 1]] <- temp[k, 2]
  }

  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  book <- data.frame(plot.num = c(t(plan.id)), row, col,
                     geno = c(t(plan)), stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]
  
  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1)
    book <- book[sort(book$plot.num, index.return = TRUE)$ix, ]
  
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
  
}
