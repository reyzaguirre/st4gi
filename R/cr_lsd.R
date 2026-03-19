#' Latin Square Design
#'
#' This function creates the fieldbook and fieldplan for a LSD.
#' @param geno The list of genotypes.
#' @param fillby Allocate the plots by \code{"rows"} or \code{"columns"},
#' default \code{"rows"}.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' cr.lsd(c("A", "B", "C"))
#' cr.lsd(c("A", "B", "C", "D", "E"))
#' @export

cr.lsd <- function(geno,
                   fillby = c('rows', 'columns'),
                   serpentine = c("yes", "no")) {
  
  # Match arguments
  
  fillby <- match.arg(fillby)
  serpentine <- match.arg(serpentine)
  
  # Error messages
  
  ng <- length(geno)
  
  if (ng < 2)
    stop("Include at least 2 genotypes.")
  
  # Fieldplan array
  
  plan <- array(dim = c(ng, ng))
  
  # Create the latin square
  
  plan[1, ] <- 1:ng

  for (i in 2:ng)
    plan[i, ] <- (plan[i - 1, ] + 1) %% ng
  
  plan[plan == 0] <- ng
  
  # Randomize rows and columns
  
  plan <- plan[sample(1:ng), ]
  plan <- plan[, sample(1:ng)]
  
  # Include genotypes at random
  
  plan <- array(geno[plan], dim = c(ng, ng))
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:ng)
  colnames(plan) <- paste("col", 1:ng)

  # Create fielbook

  row <- as.integer(gl(ng, ng))
  col <- rep(1:ng, ng)

  plan.id <- fp(ng, ng, fillby, serpentine)
  
  plot <- c(t(plan.id))
  
  book <- data.frame(plot, row, col, geno = c(t(plan)),
                     stringsAsFactors = FALSE)

  # Sort by plot number
  
  book <- book[order(book$plot), ]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
  
}
