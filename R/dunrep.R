#' Unreplicated design
#'
#' This function creates the fieldbook and fieldplan for an unreplicated design
#' with genotypes randomly allocated on a field.
#' @param geno The list of genotypes.
#' @param nr Number of rows.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' dunrep(1:100, 5)
#' dunrep(1:100, 7)
#' @export

dunrep <- function(geno, nr) {
  
  # Dimensions
  
  ng <- length(geno)      # Number of genotypes
  nc <- ceiling(ng / nr)  # Number of columns
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))

  # Include genotypes at random
  
  geno <- sample(geno)

  k <- 1

  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- geno[k]
      k <- k + 1
    }

  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Create fielbook
  
  row <- as.integer(gl(dim(plan)[1], dim(plan)[2]))
  col <- rep(1:dim(plan)[2], dim(plan)[1])
  
  book <- data.frame(plot = 1:(dim(plan)[1] * dim(plan)[2]),
                     row, col, geno = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]

  # Return
  
  list(plan = plan, book = book)
}
