#' ANOVA for MET with a RCBD
#'
#' Fit an analysis of variance model for a multi environment trial (MET) with a RCBD
#' in each environment.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default. Genotypes and environments are considered as fixed
#' factors while the blocks are considered as random and nested into the environments.
#' @return It returns the ANOVA table.
#' @author Raul Eyzaguirre.
#' @examples
#' aov.met("y", "geno", "env", "rep", met8x12)
#' @importFrom stats anova
#' @export

aov.met <- function(trait, geno, env, rep, dfr, maxp = 0.1) {

  # Everything as character

  dfr[, geno] <- as.character(dfr[, geno])
  dfr[, env] <- as.character(dfr[, env])
  dfr[, rep] <- as.character(dfr[, rep])

  # Check data
  
  lc <- ck.f(trait, c(geno, env), rep, dfr)

  # Estimate missing values and report errors from mve.met
  
  trait.est <- paste0(trait, ".est")

  if (lc$nt.0 > 0 | lc$nrep == 1 | lc$nt.mult > 0 | lc$nmis > 0 |
      lc$nmis.fac > 0 | lc$nl[1] < 2 | lc$nl[2] < 2) {
    dfr[, trait] <- mve.met(trait, geno, env, rep, dfr, maxp, tol = 1e-06)[, trait.est]
    warning(paste0("The data set is unbalanced, ",
                   format(lc$pmis * 100, digits = 3),
                   "% missing values estimated."))
  }

  # ANOVA

  model <- aov(dfr[, trait] ~ dfr[, geno] + dfr[, env]
               + dfr[, rep] %in% dfr[, env] + dfr[, geno]:dfr[, env])
  model$terms[[2]] <- trait
  
  at <- anova(model)
  
  rownames(at)[1:4] <- c(geno, env, paste0(rep, "(", env, ")"), paste0(geno, ":", env))
  
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
