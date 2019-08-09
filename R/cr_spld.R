#' Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Split-Plot 
#' or a Split-Split-Plot design with a RCBD for the whole plots.
#' @param fnames Factors' names for two or three factors. See details.
#' @param flevels A list with the factors' levels.
#' @param nb Number of blocks for the first factor.
#' @param nc Number of columns. See details.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @details Two or three factors must be included. The first factor goes to the
#' whole plots, the second to the sub-plots and the third, if any, to the sub-sub-plots.
#' The number of columns \code{"nc"} is the maximum available number of columns on
#' the field for the sub-plots or the sub-sub-plots if a third factor is provided.
#' Default is the number of levels of the last factor.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' A <- paste0("a", 1:4)
#' B <- paste0("b", 1:3)
#' cr.spld(c("A", "B"), list(A, B), 3)
#' @export

cr.spld <- function(fnames, flevels, nb, nc = NULL, serpentine = c("yes", "no")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)
  
  # Number of factors
  
  nf <- length(fnames)
  
  # Number of levels
  
  nl <- sapply(flevels, length)
  
  # Error messages
  
  if (nf < 2 | nf > 3)
    stop("Include 2 or 3 factors.")
  
  if (nf != length(flevels))
    stop("Number of factors' names does not match with the list of factors' levels")
  
  for (i in 1:nf)
    if (nl[i] < 2)
      stop(paste("Include at least 2 levels for factor", i))
  
  if (nb < 2)
    stop("Include at least 2 blocks.")
  
  # Number of rows for each plot or sub-plot
  
  if (is.null(nc))
    nc <- gnc(nl[nf])
  
  nr <- ceiling(nl[nf] / nc)

  # Fieldplan array
  
  plan.id <- fp(nr, nc, serpentine)
  
  if (nf == 2)
    plan <- array(dim = c(nr, nc, nl[1], nb))
  if (nf == 3)
    plan <- array(dim = c(nr, nc, nl[2], nl[1], nb))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
  if (nf == 2) {
    dimnames(plan)[[3]] <- paste("plot", 1:nl[1])
    dimnames(plan)[[4]] <- paste("block", 1:nb)
  }
  if (nf == 3) {
    dimnames(plan)[[3]] <- paste("subplot", 1:nl[2])
    dimnames(plan)[[4]] <- paste("plot", 1:nl[1])
    dimnames(plan)[[5]] <- paste("block", 1:nb)
  }
  
  # Include treatments at random

  for (m in 1:nb) {
    ranf1 <- sample(flevels[[1]])
    for (l in 1:nl[1]) {
      ranf2 <- sample(flevels[[2]])
      ranf12 <- paste(ranf1[l], ranf2, sep = ":-p")
      if (nf == 2)
        plan[, , l, m] <- array(ranf12[plan.id], c(nr, nc))
      if (nf == 3)
        for (k in 1:nl[2]) {
          ranf3 <- sample(flevels[[3]])
          ranf123 <- paste(ranf12[k], ranf3, sep = ":-p")
          plan[, , k, l, m] <- array(ranf123[plan.id], c(nr, nc))
        }
    }
  }
  
  # Number of smallest plot that is divided (np) and
  # number of smallest sub plots inside a plot (nsp)
  
  np <- nb * nl[1]
  nsp <- nr * nc
  if (nf == 3) {
    np <- np * nl[2]
    nsp <- nsp * nl[2]
  }
  
  # Create fielbook
  
  row <- rep(as.integer(gl(nr, nc)), np)
  col <- rep(1:nc, nr * np)
  plot.num <- as.integer(gl(nl[1] * nb, nsp))
  block <- as.integer(gl(nb, nl[1] * nsp))
  
  treat <- NULL
  subplot.num <- NULL
  subsubplot.num <- NULL
  
  if (nf == 3)
    subplot.num <- as.integer(gl(np, nr * nc))

  cont <- 0
  for (m in 1:nb)
    for (l in 1:nl[1]) {
      if (nf == 2) {
        treat <- c(treat, c(t(plan[, , l, m])))
        subplot.num <- c(subplot.num, c(t(plan.id)) + nl[2] * cont)
        cont <- cont + 1
      }
      if (nf == 3)
        for (k in 1:nl[2]) {
          treat <- c(treat, c(t(plan[, , k, l, m])))
          subsubplot.num <- c(subsubplot.num, c(t(plan.id)) + nl[3] * cont)
          cont <- cont + 1
        }
    }

  if (nf == 2)
    book <- data.frame(block, plot.num, subplot.num, row, col, treat,
                       stringsAsFactors = FALSE)
  if (nf == 3)
    book <- data.frame(block, plot.num, subplot.num, subsubplot.num, row, col, treat,
                       stringsAsFactors = FALSE)
  book <- book[!is.na(book$treat), ]
  
  # Add columns for factor levels
  
  book$f1 <- unlist(c(data.frame(sapply(book$treat, strsplit, ":-p"),
                                 stringsAsFactors = FALSE)[1, ]))

  book$f2 <- unlist(c(data.frame(sapply(book$treat, strsplit, ":-p"),
                                 stringsAsFactors = FALSE)[2, ]))

  if (nf == 3)
    book$f3 <- unlist(c(data.frame(sapply(book$treat, strsplit, ":-p"),
                                   stringsAsFactors = FALSE)[3, ]))

  nc.book <- length(colnames(book))
  colnames(book)[(nc.book - nf + 1):nc.book] <- fnames

  # Replace characters for treatment names
  
  plan <- gsub(":-p", "_", plan)
  book$treat <- gsub(":-p", "_", book$treat)

  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1 & nf == 2)
    book <- book[sort(book$subplot.num, index.return = TRUE)$ix, ]
  if (serpentine == 'yes' & nr > 1 & nf == 3)
    book <- book[sort(book$subsubplot.num, index.return = TRUE)$ix, ]
  
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
  
}
