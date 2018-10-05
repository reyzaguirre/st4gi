#' Alpha (0,1) design
#'
#' This function creates the fieldbook and fieldplan for an alpha (0,1) design.
#' @param geno The list of genotypes.
#' @param nrep Number of replications.
#' @param k Block size.
#' @param nc Number of available columns on the field.
#' @details The genotypes are randomly allocated on a field following an alpha
#' (0,1) design. In this design each replication is a complete block and is
#' divided into \code{s} incomplete blocks of size \code{k}. For any pair of
#' genotypes, the number of times they concur into the same incomplete block
#' is 0 or 1 (hence the name of the design). There are 4 options for this design:
#' \itemize{
#'  \item \code{nrep = 2} and \code{k <= s}.
#'  \item \code{nrep = 3}, \code{s} odd, and \code{k <= s}.
#'  \item \code{nrep = 3}, \code{s} even, and \code{k <= s-1}.
#'  \item \code{nrep = 4}, \code{s} even and not a multiple of 3, and \code{k <= s}.
#' }
#' The replications and blocks inside replications are disposed alongside the rows.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' cr.a01d(1:100, 2, 5, 40)
#' cr.a01d(1:100, 3, 5, 28)
#' @export

cr.a01d <- function(geno, nrep, k, nc) {
  
  # Number of genotypes and blocks
  
  ng <- length(geno)
  s <- ng / k
  
  if (ng %% k != 0)
    stop("The size of the blocks is not appropriate. The number of genotypes must be a multiple of k.")

  # Randomize list of genotypes

  geno <- geno[sample(1:ng)]
  
  # Actual number of columns
  
  nc <- floor(nc / k) * k
  
  # Design generators
  
  if (nrep == 2)
    if (k <= s) {
      alpha <- matrix(0, nrow = k, ncol = nrep)
      alpha[, 2] <- 0:(k - 1)
    } else {
      stop("With 2 replications k must be less than or equal to the number of blocks.")
    }
  
  if (nrep == 3 & s %% 2 != 0)
    if (k <= s) {
      alpha <- matrix(0, nrow = k, ncol = nrep)
      alpha[, 2] <- 0:(k - 1)
      alpha[2:k, 3] <- (s - 1):(s - k + 1)
    } else {
      stop("With 3 replications and an odd number of blocks, k must be less than or equal to the number of blocks.")
    }
  
  if (nrep == 3 & s %% 2 == 0)
    if (k < s) {
      alpha <- matrix(0, nrow = k, ncol = nrep)
      alpha[, 2] <- 0:(k - 1)
      alpha[2, 3] <- s / 2
      for (i in 3:k)
        alpha[i, 3] <- alpha[i - 2, 3] + 1
    } else {
      stop("With 3 replications and an even number of blocks, k must be less than the number of blocks.")
    }
  
  if (nrep == 4)
    if (s %% 2 != 0 & s %% 3 != 0 & k <= s) {
      alpha <- matrix(0, nrow = k, ncol = nrep)
      alpha[, 2] <- 0:(k - 1)
      alpha[2:k, 3] <- (s - 1):(s - k + 1)
      alpha[2, 4] <- (s + 1) / 2
      for (i in 3:k)
        alpha[i, 4] <- alpha[i - 2, 4] + 1
    } else {
      if (s %% 2 == 0)
        stop("With 4 replications the number of blocks must be odd.")
      if (s %% 3 == 0)
        stop("With 4 replications the number of blocks cannot be a multiple of 3.")
      if (k > s)
        stop("With 4 replications k must be less than or equal to the number of blocks.")
    }
  
  if (nrep > 4)
    stop("The maximum number of replications for this design is 4.")
  
  # Repeat desin generator s times
  
  ad <- NULL
  
  for (i in 1:nrep)
    ad <- c(ad, rep(alpha[, i], s))

  dim(ad) <- c(k, s, nrep)
  
  # Cyclic substitution

  for (i in 1:nrep)
    for (j in 2:s)
      ad[, j, i] <- (ad[, j - 1, i] + 1) %% s

  # Add s, 2s, 3s, ... to each row
  
  for (i in 2:k)
    ad[i, , ] <- ad[i, , ] + (i - 1) * s
  
  # Add 1 to get genotype numbers from 1 to ng
  
  ad <- ad + 1
  
  # Randomize genotypes inside each block
  
  for (i in 1:nrep)
    for (j in 1:s)
      ad[, j, i] <- sample(ad[, j, i])
    
  # Randomize blocks inside each replication
  
  for (i in 1:nrep)
    ad[, , i] <- ad[, sample(1:s), i]

  # Number of rows for each rep
  
  nr <- ceiling(ng / nc)
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc, nrep))

  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("rep", 1:nrep)
  
  # Allocate genotypes in the fieldplan
  
  for (i in 1:nrep) {
    sg <- c(ad[, , i])
    m <- 1
    for (j in 1:nr)
      for (l in 1:nc) {
        plan[j, l, i] <- geno[sg[m]]
        m <- m + 1
      }
  }
  
  # Create fielbook

  r <- as.integer(gl(nrep, nr * nc))
  row <- rep(as.integer(gl(nr, nc)), nrep)
  col <- rep(rep(1:nc, nr), nrep)
  block <- as.integer(gl(s, k))
  length(block) <- nr * nc
  block <- rep(block, nrep)
  
  geno <- NULL
  
  for (i in 1:nrep)
    geno <- c(geno, c(t(plan[, , i])))
  
  book <- data.frame(r, block, row, col, geno, stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  book$plot <- 1:dim(book)[1]
  book <- book[, c(6, 1, 2, 3, 4, 5)]
  rownames(book) <- 1:dim(book)[1]

  # Return
  
  list(plan = plan, book = book)
  
}
