#' Create design for a factorial experiment
#'
#' This function creates the fieldbook and fieldplan for a factorial experiment
#' following a CRD or a RCBD.
#' @param fnames Factors' names.
#' @param flevels A list with the factors' levels.
#' @param design The design, \code{crd} or \code{rcbd}.
#' @param nrep Number of replications or blocks.
#' @param nc Number of available columns on the field.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @param alongside \code{"no"} for independent blocks, or \code{"rows"}
#' or \code{"columns"} if blocks are together alongside rows or columns.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' A <- paste0("a", 1:5)
#' B <- paste0("b", 1:3)
#' C <- paste0("c", 1:2)
#' cr.f(c("A", "B", "C"), list(A, B, C), "rcbd", 3, 12)
#' @export

cr.f <- function(fnames, flevels, design = c("crd", "rcbd"), nrep, nc = NULL,
                 serpentine = c("yes", "no"), alongside = c("no", "rows", "columns")) {
  
  # Match arguments
  
  design <- match.arg(design)
  serpentine <- match.arg(serpentine)
  alongside <- match.arg(alongside)
  
  # Number of factors
  
  nf <- length(fnames)

  # Number of levels
  
  nl <- sapply(flevels, length)
  
  # Error messages
  
  if (nf < 2)
    stop("Include at least 2 factors.")
  
  if (nf != length(flevels))
    stop("Number of factors' names does not match with the list of factors' levels")
    
  for (i in 1:nf)
    if (nl[i] < 2)
      stop(paste("Include at least 2 levels for factor", i))

  # Number of treatments
  
  nt <- prod(sapply(flevels, length))

  # Treatment labels
  
  trt <- c(outer(flevels[[1]], flevels[[2]], paste, sep = ":-p"))
  
  if (nf > 2)
    for (i in 3:nf)
      trt <- c(outer(trt, flevels[[i]], paste, sep = ":-p"))
  
  # Create fielbook and fieldplan
  
  if (design == "crd")
    output <- cr.crd(trt, nrep, nc, serpentine)
  
  if (design == "rcbd")
    output <- cr.rcbd(trt, nrep, nc, serpentine, alongside)
  
  # Add columns for factor levels
  
  temp <- unlist(strsplit(output$book$geno, ":-p"))
  
  fm <- matrix(temp, nrow = length(temp) / nf, ncol = nf, byrow = TRUE)
  
  for (i in 1:nf)
    output$book[, fnames[i]] <- fm[, i]
  
  # Rename geno by treat
  
  if (design == "crd")
    colnames(output$book)[4] <- "treat"

  if (design == "rcbd")
    colnames(output$book)[5] <- "treat"

  # Replace characters for treatment names
  
  output$book$treat <- gsub(":-p", "_", output$book$treat)
  output$plan <- gsub(":-p", "_", output$plan)
    
  # Return
  
  list(plan = output$plan, book = output$book)

}
