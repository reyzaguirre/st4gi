#' ANOVA for a 2-factor factorial
#'
#' Fit an analysis of variance model for a 2-factor factorial with a CRD or RCBD
#' @param trait The trait to analyze.
#' @param A Factor A.
#' @param B Factor B.
#' @param rep The replications or blocks.
#' @param design The statistical design, \code{crd} or \code{rcbd}.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @author Raul Eyzaguirre.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default.
#' @return It returns the ANOVA table.
#' @examples
#' aov.2f("asc.dw", "geno", "treat", "rep", "crd", asc)
#' @importFrom stats anova
#' @export

aov.2f <- function(trait, A, B, rep, design = c("crd", "rcbd"), data, maxp = 0.1) {

  # match arguments
  
  design <- match.arg(design)

    # Everything as factor

  data[, A] <- factor(data[, A])
  data[, B] <- factor(data[, B])
  data[, rep] <- factor(data[, rep])

  # Check data and estimate missing values

  lc <- check.2f(trait, A, B, rep, data)

  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0) {
    data[, trait] <- mve.2f(trait, A, B, rep, design, data, maxp, tol = 1e-06)[, 5]
    warning(paste("The data set is unbalanced, ",
                  format(lc$pmis * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }

  # Error messages

  if (lc$na < 2 | lc$nb < 2)
    stop("This is not a 2-factor factorial.")

  # ANOVA
  
  if (design == "crd") {
    model <- aov(data[, trait] ~ data[, A] * data[, B])
    i <- 3
    rn <- c(A, B, paste(A, ":", B, sep = ""))
  }
  
  if (design == "rcbd") {
    model <- aov(data[, trait] ~ data[, A] * data[, B] + data[, rep])
    i <- 4
    rn <- c(A, B, rep, paste(A, ":", B, sep = ""))
  }
  
  model$terms[[2]] <- trait
  
  at <- anova(model)
  
  rownames(at)[1:i] <- rn
  
  # Correction for missing values
  
  if (lc$nmis > 0) {
    at[i + 1, 1] <- at[i + 1, 1] - lc$nmis
    at[i + 1, 3] <- at[i + 1, 2] / at[i + 1, 1]
    at[1:i, 4] <- at[1:i, 3] / at[i + 1, 3]
    at[1:i, 5] <- pf(at[1:i, 4], at[1:i, 1], at[i + 1, 1], lower.tail = FALSE)
  }
  
  # Return
  
  at
}
