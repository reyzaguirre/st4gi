#' ANOVA for MET with a RCBD
#'
#' Fit an analysis of variance model for a multi environment trial (MET) with a RCBD
#' in each environment.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @author Raul Eyzaguirre.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default. Genotypes and environments are considered as fixed
#' factors while the blocks are considered as random and nested into the environments.
#' @return It returns the ANOVA table.
#' @examples
#' aovmet("y", "geno", "env", "rep", met8x12)
#' @importFrom stats anova
#' @export

aovmet <- function(trait, geno, env, rep, data, maxp = 0.1) {

  # Everything as factor

  data[, geno] <- factor(data[, geno])
  data[, env] <- factor(data[, env])
  data[, rep] <- factor(data[, rep])

  # Check data and estimate missing values

  lc <- check.met(trait, geno, env, rep, data)

  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0) {
    data[, trait] <- mvemet(trait, geno, env, rep, data, maxp, tol = 1e-06)[, 5]
    warning(paste("The data set is unbalanced, ",
                  format(lc$pmis * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }

  # Error messages

  if (lc$ng < 2 | lc$ne < 2)
    stop("This is not a MET experiment.")

  # ANOVA

  model <- aov(data[, trait] ~ data[, geno] + data[, env]
               + data[, rep] %in% data[, env] + data[, geno]:data[, env])
  model$terms[[2]] <- trait
  
  at <- anova(model)
  
  rownames(at)[1:4] <- c(geno, env, paste(rep, "(", env, ")", sep = ""),
                         paste(geno, ":", env, sep = ""))
  
  # Correction for missing values
  
  if (lc$nmis > 0) {
    at[5, 1] <- at[5, 1] - lc$nmis
    at[5, 3] <- at[5, 2] / at[5, 1]
    at[c(1, 3, 4), 4] <- at[c(1, 3, 4), 3] / at[5, 3]
    at[c(1, 3, 4), 5] <- pf(at[c(1, 3, 4), 4], at[c(1, 3, 4), 1], at[5, 1], lower.tail = FALSE)
  }
  
  # Correction for blocks nested into environments
  
  at[2, 4] <- at[2, 3] / at[3, 3]
  at[2, 5] <- pf(at[2, 4], at[2, 1], at[3, 1], lower.tail = FALSE)  

  # Return
  
  at
}
