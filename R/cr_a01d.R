#' Alpha (0,1) design
#'
#' This function creates the fieldbook and fieldplan for an alpha (0,1) design.
#' @param geno The list of genotypes.
#' @param check The list of checks to repeat in each incomplete block (optional).
#' @param nb Number of complete blocks.
#' @param k Size for the incomplete blocks.
#' @param nc Number of available columns on the field for each complete block .
#' @param breakib If incomplete blocks should be broken in more than one
#' row, default \code{"no"}.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @param alongside \code{"no"} for independent blocks, or \code{"rows"}
#' or \code{"columns"} if blocks are together alongside rows or columns.
#' @details The genotypes are randomly allocated on a field following an alpha
#' (0,1) design. In this design each block is a complete replication that is
#' divided into \code{s} incomplete blocks of size \code{k}. For any pair of
#' genotypes, the number of times they concur into the same incomplete block
#' is 0 or 1 (hence the name of the design). There are 4 options for this design:
#' \itemize{
#'  \item \code{nb = 2} and \code{k <= s}.
#'  \item \code{nb = 3}, \code{s} odd, and \code{k <= s}.
#'  \item \code{nb = 3}, \code{s} even, and \code{k <= s-1}.
#'  \item \code{nb = 4}, \code{s} even and not a multiple of 3, and \code{k <= s}.
#' }
#' The blocks and incomplete blocks are disposed alongside the rows.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' cr.a01d(1:100, NULL, 2, 5)
#' cr.a01d(1:100, NULL, 3, 5, 28)
#' @export

cr.a01d <- function(geno, check, nb, k, nc = NULL, breakib = c('no', 'yes'),
                    serpentine = c("yes", "no"), alongside = c("no", "rows", "columns")) {
  
  # Match arguments
  
  breakib <- match.arg(breakib)
  serpentine <- match.arg(serpentine)
  alongside <- match.arg(alongside)
  
  # Number of genotypes, checks, and incomplete blocks
  
  nck <- length(check)
  ng <- length(geno)
  s <- ng / k
  
  if (ng %% k != 0)
    stop("The size of the incomplete blocks is not appropriate.
         The number of genotypes must be a multiple of k.")
  
  # Number of genotypes plus checks
  
  ngc <- ng + nck * s

  # Randomize list of genotypes

  geno <- geno[sample(1:ng)]
  
  # Actual number of columns
  
  if (is.null(nc))
    nc <- round(sqrt(ngc))
  
  if (floor(nc / (k + nck)) > 0 & breakib == 'no')
    nc <- floor(nc / (k + nck)) * (k + nck)

  # Design generators
  
  if (nb == 2)
    if (k <= s) {
      alpha <- matrix(0, nrow = k, ncol = nb)
      alpha[, 2] <- 0:(k - 1)
    } else {
      stop("With 2 complete blocks, k must be less than or equal
           to the number of incomplete blocks.")
    }
  
  if (nb == 3 & s %% 2 != 0)
    if (k <= s) {
      alpha <- matrix(0, nrow = k, ncol = nb)
      alpha[, 2] <- 0:(k - 1)
      alpha[2:k, 3] <- (s - 1):(s - k + 1)
    } else {
      stop("With 3 complete blocks and an odd number of incomplete blocks,
           k must be less than or equal to the number of incomplete blocks.")
    }
  
  if (nb == 3 & s %% 2 == 0)
    if (k < s) {
      alpha <- matrix(0, nrow = k, ncol = nb)
      alpha[, 2] <- 0:(k - 1)
      alpha[2, 3] <- s / 2
      for (i in 3:k)
        alpha[i, 3] <- alpha[i - 2, 3] + 1
    } else {
      stop("With 3 complete blocks and an even number of incomplete blocks,
           k must be less than the number of incomplete blocks.")
    }
  
  if (nb == 4)
    if (s %% 2 != 0 & s %% 3 != 0 & k <= s) {
      alpha <- matrix(0, nrow = k, ncol = nb)
      alpha[, 2] <- 0:(k - 1)
      alpha[2:k, 3] <- (s - 1):(s - k + 1)
      alpha[2, 4] <- (s + 1) / 2
      for (i in 3:k)
        alpha[i, 4] <- alpha[i - 2, 4] + 1
    } else {
      if (s %% 2 == 0)
        stop("With 4 complete blocks the number of incomplete blocks must be odd.")
      if (s %% 3 == 0)
        stop("With 4 complete blocks the number of incomplete blocks
             cannot be a multiple of 3.")
      if (k > s)
        stop("With 4 complete blocks k must be less than or equal to
             the number of incomplete blocks.")
    }
  
  if (nb > 4)
    stop("The maximum number of complete blocks for this design is 4.")
  
  # Repeat desin generator s times
  
  ad <- NULL
  
  for (i in 1:nb)
    ad <- c(ad, rep(alpha[, i], s))

  dim(ad) <- c(k, s, nb)
  
  # Cyclic substitution

  for (i in 1:nb)
    for (j in 2:s)
      ad[, j, i] <- (ad[, j - 1, i] + 1) %% s

  # Add s, 2s, 3s, ... to each row
  
  for (i in 2:k)
    ad[i, , ] <- ad[i, , ] + (i - 1) * s
  
  # Add 1 to get genotype numbers from 1 to ng
  
  ad <- ad + 1
  
  # Add checks
  
  if (nck > 0) {
    ad2 <- array(NA, c(k + nck, s, nb))
    ad2[1:k, 1:s, 1:nb] <- ad
    for (i in 1:nck)
      ad2[k + i, , ] <- ng + i
    ad <- ad2
  }
   
  # Randomize genotypes inside each incomplete block
  
  for (i in 1:nb)
    for (j in 1:s)
      ad[, j, i] <- sample(ad[, j, i])
    
  # Randomize incomplete blocks inside each complete block
  
  for (i in 1:nb)
    ad[, , i] <- ad[, sample(1:s), i]

  # Number of rows for each complete block
  
  nr <- ceiling(ngc / nc)
  
  # Fieldplan array
  
  plan.id <- t(array(1:(nr*nc), dim = c(nc, nr)))
  plan.id.sb <- t(array(c(sapply(1:s, rep, k + nck)), dim = c(nc, nr)))
  
  if (serpentine == 'yes' & nr > 1)
    for (i in seq(2, nr, 2)) {
      plan.id[i, ] <- sort(plan.id[i, ], decreasing = TRUE)
      plan.id.sb[i, ] <- sort(plan.id.sb[i, ], decreasing = TRUE)
    }

  plan <- array(dim = c(nr, nc, nb))

  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("block", 1:nb)
  
  # Allocate genotypes in the fieldplan
  
  geno <- c(geno, check)
  
  for (l in 1:nb) {
    sg <- geno[c(ad[, , l])]
    plan[, , l] <- array(sg[plan.id], c(nr, nc))
  }
  
  # Create fielbook

  block <- as.integer(gl(nb, nr * nc))
  row <- rep(as.integer(gl(nr, nc)), nb)
  col <- rep(rep(1:nc, nr), nb)

  geno <- NULL
  plot.num <- NULL
  sub.block <- NULL
  
  for (i in 1:nb) {
    geno <- c(geno, c(t(plan[, , i])))
    plot.num <- c(plot.num, c(t(plan.id)) + ngc * (i - 1))
    sub.block <- c(sub.block, c(t(plan.id.sb)))
  }
  
  book <- data.frame(plot.num, block, sub.block, row, col, geno, stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]

  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1)
    book <- book[sort(book$plot.num, index.return = TRUE)$ix, ]
  
  rownames(book) <- 1:dim(book)[1]

  # Change row and column numbers if required
  
  if (alongside == "rows") {
    plan <- t(apply(plan, 1, rbind))
    colnames(plan) <- paste("col", 1:dim(plan)[2])
    book$col <- book$col + (book$block - 1) * nc
  }
  
  if (alongside == "columns") {
    plan <- apply(plan, 2, rbind)
    rownames(plan) <- paste("row", 1:dim(plan)[1])
    book$row <- book$row + (book$block - 1) * nr
  }  

  # Return
  
  list(plan = plan, book = book)
  
}
