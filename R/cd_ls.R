#' Latin Square Design
#'
#' This function creates the fieldbook and fieldplan for a LSD.
#' @param geno The list of genotypes.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field following a LSD.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' cd.ls(c("A", "B", "C"))
#' cd.ls(c("A", "B", "C", "D", "E"))
#' @export

cd.ls <- function(geno) {
  
  # Error messages
  
  ng <- length(geno)
  if (ng < 2)
    stop("Include at least 2 genotypes.")
  
  # Fieldplan array
  
  plan <- array(dim = c(ng, ng))
  plan[1, ] <- 1:ng
  for (i in 2:ng)
    plan[i, ] <- plan[i - 1, ] + 1
  plan <- plan %% ng
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
  book <- data.frame(plot = 1:(ng * ng), row, col,
                     geno = c(t(plan)), stringsAsFactors = F)

  # Return
  
  list(plan = plan, book = book)
}
