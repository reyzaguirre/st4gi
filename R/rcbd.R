#' ANOVA for a RCBD
#'
#' Fit an analysis of variance model for a RCBD.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param block The blocks.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param anova Logical, if TRUE the ANOVA table is shown.
#' @author Raul Eyzaguirre
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default.
#' @return If \code{anova} is {TRUE} it returns and shows the ANOVA table.
#' If \code{anova} is {FALSE} it returns the ANOVA table and some other components as the
#' estimated missing values, but nothing is shown.
#' @examples
#' # The data
#' head(pjpz09)
#' str(pjpz09)
#'
#' # Get a copy with some missing values for trw
#' temp <- pjpz09
#' temp[c(10, 20, 30), "trw"] <- NA
#' 
#' # Run ANOVA for trw
#' rcbd("trw", "geno", "rep", temp)
#' @export

rcbd <- function(trait, treat, block, data, maxp = 0.1, anova = TRUE) {

  # Everything as factor

  data[, treat] <- factor(data[, treat])
  data[, block] <- factor(data[, block])

  # Check data and estimate missing values

  lc <- checkdata01(trait, treat, data)

  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0) {
    est.data <- mveb(trait, treat, block, data, maxp, tol = 1e-06)
    data[, trait] <- est.data$new.data[, 4]
    nmis <- est.data$est.num
    warning(paste("The data set is unbalanced, ",
                  format(est.data$est.prop * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  } else {
    nmis <- 0
  }

  # ANOVA

  model <- aov(data[, trait] ~ data[, treat] + data[, block])
  model$terms[[2]] <- trait
  
  at <- anova(model)
  
  rownames(at)[1:2] <- c(treat, block)
  
  # Correction for missing values
  
  at[3, 1] <- at[3, 1] - nmis
  at[3, 3] <- at[3, 2] / at[3, 1]
  at[c(1, 2), 4] <- at[c(1, 2), 3] / at[3, 3]
  at[c(1, 2), 5] <- pf(at[c(1, 2), 4], at[c(1, 2), 1], at[3, 1], lower.tail = FALSE)
  
  # Return

  output <- list (at = at, lc = lc, est.data = est.data)
  
  if(anova == TRUE) output$at else invisible(output)
}
