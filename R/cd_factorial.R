#' Create design for a factorial experiment
#'
#' This function creates the fieldbook and fieldplan for a factorial experiment
#' with 2 to 5 factors following a CRD or a RCBD.
#' @param A The levels of factor A.
#' @param B The levels of factor B.
#' @param C The levels of factor C.
#' @param D The levels of factor D.
#' @param E The levels of factor E.
#' @param design The design, a crd or a rcbd.
#' @param nrep Number of replications or blocks.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The treatments are randomly allocated on a field following a CRD or a RCBD.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' A <- paste("a", 1:5, sep = "")
#' B <- paste("b", 1:3, sep = "")
#' cd.factorial(A, B, design = "rcbd", nrep = 3, nc = 12)
#' @export

cd.factorial <- function(A, B, C = NULL, D = NULL, E = NULL,
                         design = c("crd", "rcbd"), nrep, nc) {
  
  # Match arguments
  
  design <- match.arg(design)
  
  # Number of factors
  
  nf <- sum(!sapply(list(A, B, C, D, E), is.null))
  
  # Factor levels as characters
  
  A <- as.character(A)
  B <- as.character(B)
  C <- as.character(C)
  D <- as.character(D)
  E <- as.character(E)
  
  # Error messages
  
  nla <- length(A)
  nlb <- length(B)
  nlc <- length(C)
  nld <- length(D)
  nle <- length(E)
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  if (nla < 2)
    stop("Include at least 2 levels for factor A.")

  if (nlb < 2)
    stop("Include at least 2 levels for factor B.")

  if (nlc == 1)
    stop("Include at least 2 levels for factor C.")

  if (nld == 1)
    stop("Include at least 2 levels for factor D.")

  if (nle == 1)
    stop("Include at least 2 levels for factor E.")
  
  # Number of treatments
  
  if (nlc == 0) nlc <- 1
  if (nld == 0) nld <- 1
  if (nle == 0) nle <- 1
  
  nt <- nla * nlb * nlc * nld * nle

  # Treatment labels
  
  trt <- c(outer(A, B, paste, sep = "_"))
  
  if (nlc > 1)
    trt <- c(outer(trt, C, paste, sep = "_"))
  
  if (nld > 1)
    trt <- c(outer(trt, D, paste, sep = "_"))
  
  if (nle > 1)
    trt <- c(outer(trt, E, paste, sep = "_"))
  
  # Create fielbook and fieldplan
  
  if (design == "crd")
    output <- cd.cr(trt, nrep, nc)
  
  if (design == "rcbd")
    output <- cd.rcb(trt, nrep, nc)
  
  # Add columns to fielbook
  
  temp <- unlist(strsplit(output$book$geno, "_"))
  
  fm <- matrix(temp, nrow = length(temp) / nf, ncol = nf, byrow = TRUE)
  
  fn <- c("A", "B", "C", "D", "E")
  
  for (i in 1:nf)
    output$book[, fn[i]] <- fm[, i]
  
  if (design == "crd")
    colnames(output$book)[4] <- "treat"

  if (design == "rcbd")
    colnames(output$book)[5] <- "treat"

  # Return
  
  list(plan = output$plan, book = output$book)
}
