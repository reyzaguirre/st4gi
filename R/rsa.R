#' Regression Stability Analysis
#'
#' Function to run the regression stability analysis (Yates and Cochran, 1938,
#' Finlay and Wilkinson, 1963).
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @details The regression stability analysis is evaluated with a balanced data set.
#' If data is unbalanced, missing values are estimated up to an specified maximum proportion,
#' 10\% by default. For the ANOVA table, genotypes and environments are considered as fixed
#' factors while the blocks are considered as random and nested into the environments.
#' To run a regression stability analysis you need a set of genotypes evaluated in a set of
#' environments. At least 3 genotypes or environments are needed. In a regression stability
#' analysis for genotypes grown at several environments, for each genotype a simple linear
#' regression of individual yield (Y) on the mean yield of all genotypes for each environment
#' (X) is fitted. In a similar way, for each environment a simple linear regression of
#' individual yield (Y) on the mean yield of all environments for each genotype (X) is fitted.
#' In both cases the X values are centered on zero, so the intercepts of the models correspond to
#' the means of the genotypes or environments.
#' @return It returns the regression stability analysis decomposition of the GxE interaction
#' for genotypes and environments (Heterogeneity among regressions and deviation from regression),
#' the coefficient of variation, and the following regression stability measures for genotypes
#' and environments:
#' \itemize{
#' \item \code{a} the intercept. 
#' \item \code{b} the slope.
#' \item \code{se} the standard error for the slope.
#' \item \code{MSe} the mean square error.
#' \item \code{MSentry} the variance of the genotype means across environments and the environment
#' means across genotypes.
#' \item \code{MSinter} the variance of the genotype interaction effects across environments and the
#' environment interaction effects across genotypes.
#' }
#' @author Raul Eyzaguirre.
#' @references
#' Finlay, K. W., and Wilkinson, G. N. (1963). The Analysis of Adaption in a Plant-Breeding Programme.
#' Aust. J. Agric. Res. 14: 742-754.
#'
#' Yates, F., and Cochran, W. G. (1938). The Analysis of Group Experiments.
#' J. Agric. Sci. 28: 556-580.
#' @examples
#' rsa("y", "geno", "env", "rep", met8x12)
#' @importFrom stats coef lm summary.lm
#' @export

rsa <- function(trait, geno, env, rep, dfr, maxp = 0.1) {
  
  # Error messages
  
  lc <- ck.f(trait, c(geno, env), rep, dfr)
  
  if (lc$nl[1] < 3 & lc$nl[2] < 3)
    stop("You need at least 3 genotypes or 3 environments for regression stability analysis.")

  # Compute ANOVA and report errors from mve.met
   
  at <- suppressWarnings(aov.met(trait, geno, env, rep, dfr, maxp))

  # Estimate missing values

  trait.est <- paste0(trait, ".est")
  
  if (lc$nt.0 > 0 | lc$nrep == 1 | lc$nt.mult > 0 | lc$nmis > 0 | lc$nmis.fac > 0) {
    dfr[, trait] <- mve.met(trait, geno, env, rep, dfr, maxp, tol = 1e-06)[, trait.est]
    warning(paste0("The data set is unbalanced, ",
                   format(lc$pmis * 100, digits = 3),
                   "% missing values estimated."))
  }
  
  # Some statistics

  int.mean <- tapply(dfr[, trait], list(dfr[, geno], dfr[, env]), mean, na.rm = TRUE)
  overall.mean <- mean(int.mean, na.rm = TRUE)
  env.mean <- apply(int.mean, 2, mean, na.rm = TRUE)
  geno.mean <- apply(int.mean, 1, mean, na.rm = TRUE)
  
  # Regression-stability for genotypes

  a <- NULL        # linear regression intercept
  b <- NULL        # linear regression slope
  se <- NULL       # slope standard error
  MSe <- NULL      # error mean square
  MSentry <- NULL  # variance of the means
  MSinter <- NULL  # variance of the interaction effects
  ssr <- NULL      # residual sum of squares

  for (i in 1:lc$nl[1]) {
    modelo <- lm(int.mean[i, ] ~ I(env.mean - overall.mean))
    a[i] <- coef(modelo)[1]
    b[i] <- coef(modelo)[2]
    se[i] <- summary.lm(modelo)$coefficients[2, 2]
    MSe[i] <- anova(modelo)[2, 3]
    MSentry[i] <- sum((int.mean[i, ] - geno.mean[i])^2) / (lc$nl[2] - 1)
    MSinter[i] <- sum((int.mean[i, ] - geno.mean[i] - env.mean + overall.mean)^2) / (lc$nl[2] - 1)
    ssr[i] <- anova(modelo)[2, 2] * lc$nrep
  }
  stab.geno <- cbind(a, b, se, MSe, MSentry, MSinter)
  row.names(stab.geno) <- row.names(int.mean)
  
  if (lc$nl[2] > 2) {
    drg.sc <- sum(ssr)
    hrg.sc <- at[4, 2] - drg.sc
    hrg.gl <- lc$nl[1] - 1
    drg.gl <- (lc$nl[1] - 1) * (lc$nl[2] - 1) - hrg.gl
    drg.cm <- drg.sc / drg.gl
    hrg.cm <- hrg.sc / hrg.gl
    drg.f <- drg.cm / at[5, 3]
    hrg.f <- hrg.cm / drg.cm
    drg.p <- pf(drg.f, drg.gl, at[5, 1], lower.tail = FALSE)
    hrg.p <- pf(hrg.f, hrg.gl, drg.gl, lower.tail = FALSE)
  } else {
    drg.sc <- NA
    hrg.sc <- NA
    drg.gl <- NA
    hrg.gl <- NA
    drg.cm <- NA
    hrg.cm <- NA
    drg.f <- NA
    hrg.f <- NA
    drg.p <- NA
    hrg.p <- NA
  }

  # Regression-stability for environments

  a <- NULL        # linear regression intercept
  b <- NULL        # linear regression slope
  se <- NULL       # slope standard error
  MSe <- NULL      # error mean square
  MSentry <- NULL  # variance of the means
  MSinter <- NULL  # variance of the interaction effects
  ssr <- NULL      # residual sum of squares
  
  for (i in 1:lc$nl[2]) {
    modelo <- lm(int.mean[, i] ~ I(geno.mean - overall.mean))
    a[i] <- coef(modelo)[1]
    b[i] <- coef(modelo)[2]
    se[i] <- summary.lm(modelo)$coefficients[2, 2]
    MSe[i] <- anova(modelo)[2, 3]
    MSentry[i] <- sum((int.mean[, i] - env.mean[i])^2) / (lc$nl[1] - 1)
    MSinter[i] <- sum((int.mean[, i] - env.mean[i] - geno.mean + overall.mean)^2) / (lc$nl[1] - 1)
    ssr[i] <- anova(modelo)[2, 2] * lc$nrep
  }
  stab.env <- cbind(a, b, se, MSe, MSentry, MSinter)
  row.names(stab.env) <- colnames(int.mean)

  if (lc$nl[1] > 2) {
    dre.sc <- sum(ssr)
    hre.sc <- at[4, 2] - dre.sc
    hre.gl <- lc$nl[2] - 1
    dre.gl <- (lc$nl[1] - 1) * (lc$nl[2] - 1) - hre.gl
    dre.cm <- dre.sc / dre.gl
    hre.cm <- hre.sc / hre.gl
    dre.f <- dre.cm / at[5, 3]
    hre.f <- hre.cm / dre.cm
    dre.p <- pf(dre.f, dre.gl, at[5, 1], lower.tail = FALSE)
    hre.p <- pf(hre.f, hre.gl, dre.gl, lower.tail = FALSE)
  } else {
    dre.sc <- NA
    hre.sc <- NA
    dre.gl <- NA
    hre.gl <- NA
    dre.cm <- NA
    hre.cm <- NA
    dre.f <- NA
    hre.f <- NA
    dre.p <- NA
    hre.p <- NA
  }

  # ANOVA plus regression stability

  cv <- sqrt(at[5, 3]) / abs(overall.mean) * 100

  fileaux <- at[5, ]
  at[5, ] <- c(hrg.gl, hrg.sc, hrg.cm, hrg.f, hrg.p)
  at[6, ] <- c(drg.gl, drg.sc, drg.cm, drg.f, drg.p)
  at[7, ] <- c(hre.gl, hre.sc, hre.cm, hre.f, hre.p)
  at[8, ] <- c(dre.gl, dre.sc, dre.cm, dre.f, dre.p)
  at[9, ] <- fileaux
  row.names(at) <- c("G", "E", "R:E", "GxE", "- Het.Regr.G", "- Dev.Regr.G",
                     "- Het.Regr.E", "- Dev.Regr.E", "Residuals")

  # Return

  list(ANOVA = at, CV = cv, Stability_for_genotypes = stab.geno,
       Stability_for_environments = stab.env)
  
}
