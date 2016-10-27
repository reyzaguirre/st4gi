#' ANOVA for a RCBD
#'
#' Fit an analysis of variance model for a RCBD.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param rep The replications.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @author Raul Eyzaguirre.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default.
#' @return It returns ANOVA table.
#' @examples
#' ## Get a copy with some missing values for trw and run ANOVA
#' temp <- pjpz09
#' temp[c(10, 20, 30), "trw"] <- NA
#' aov.rcbd("trw", "geno", "rep", temp)
#' @export

aov.rcbd <- function(trait, treat, rep, data, maxp = 0.1) {

  # Everything as factor

  data[, treat] <- factor(data[, treat])
  data[, rep] <- factor(data[, rep])

  # Check data and estimate missing values

  lc <- check.rcbd(trait, treat, rep, data)

  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0 | lc$c4 == 0) {
    data[, trait] <- mve.rcbd(trait, treat, rep, data, maxp, tol = 1e-06)[, 4]
    warning(paste("The data set is unbalanced, ",
                  format(lc$pmis * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }

  # ANOVA

  model <- aov(data[, trait] ~ data[, treat] + data[, rep])
  model$terms[[2]] <- trait
  
  at <- anova(model)
  
  rownames(at)[1:2] <- c(treat, rep)
  
  # Correction for missing values
  
  if (lc$nmis > 0) {
    at[3, 1] <- at[3, 1] - lc$nmis
    at[3, 3] <- at[3, 2] / at[3, 1]
    at[c(1, 2), 4] <- at[c(1, 2), 3] / at[3, 3]
    at[c(1, 2), 5] <- pf(at[c(1, 2), 4], at[c(1, 2), 1], at[3, 1], lower.tail = FALSE)
  }
  
  # Return

  at
}
