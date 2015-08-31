#' Regression Stability Analysis
#'
#' Function to run the regression stability analysis (Yates and Cochran, 1938,
#' Finlay and Wilkinson, 1963). This implementation follows the formulas of
#' Eberhart and Russell (1966).
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @author Raul Eyzaguirre
#' @details The regression stability analysis is evaluated with a balanced data set.
#' If data is unbalanced, missing values are estimated up to an specified maximum proportion,
#' 10\% by default. For the ANOVA table, genotypes and environments are considered as fixed
#' factors while the blocks are considered as random and nested into the environments.
#' To run a regression stability analysis you need a set of genotypes evaluated in a set of
#' environments. At least 3 genotypes or environments are needed. In a regression stability
#' analysis for genotypes grown at several environments, for each genotype a simple linear
#' regression of individual yield (Y) on the mean yield of all genotypes for each environment
#' (X) is fitted. In a similar way, for each environment, a simple linear regression of
#' individual yield (Y) on the mean yield of all environments for each genotype (X) is fitted.
#' @return It returns the regression stability analysis decomposition of the GxE interaction
#' for genotypes and environments, the coefficient of variation, and the regression stability
#' measures for genotypes and environments. 
#' @references
#' Eberhart, S. A. and Russell, W. A. (1966). Stability Parameters for Comparing Varieties.
#' Crop Sci. 6: 36-40.
#'
#' Finlay, K. W., and Wilkinson, G. N. (1963). The Analysis of Adaption in a Plant-Breeding Programme.
#' Aust. J. Agric. Res. 14: 742-754.
#'
#' Yates, F., and Cochran, W. G. (1938). The Analysis of Group Experiments.
#' J. Agric. Sci. 28: 556-580.
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Run regression stability analysis
#' rsa("y", "geno", "env", "rep", met8x12)
#' @export

rsa <- function(trait, geno, env, rep, data, maxp = 0.1) {
  
  # Error messages
  
   geno.num <- nlevels(data[, geno])
   env.num <- nlevels(data[, env])
  
  if (geno.num == 2 & env.num == 2)
    stop("You need at least 3 genotypes or 3 environments for regression stability analysis.")

  # Compute ANOVA
  
  at <- suppressWarnings(aovmet(trait, geno, env, rep, data, maxp))

  # Check data and estimate missing values
  
  lc <- checkdata02(trait, geno, env, data)
  
  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0) {
    est.data <- mvemet(trait, geno, env, rep, data, maxp, tol = 1e-06)
    data[, trait] <- est.data$new.data[, 5]
    warning(paste("The data set is unbalanced, ",
                  format(est.data$est.prop * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }
  
  # Some statistics

  int.mean <- tapply(data[, trait], list(data[, geno], data[, env]), mean, na.rm = TRUE)
  overall.mean <- mean(int.mean, na.rm = TRUE)
  env.mean <- apply(int.mean, 2, mean, na.rm = TRUE)
  geno.mean <- apply(int.mean, 1, mean, na.rm = TRUE)
  
  # Regression-stability for genotypes

  a <- NULL        # linear regression intercept
  b <- NULL        # linear regression slope
  se <- NULL       # slope standard error
  ms_e <- NULL     # error mean square
  ms_gxe <- NULL   # variance of the interaction effects for genotype i

  for (i in 1:geno.num) {
    modelo <- lm(int.mean[i, ] ~ env.mean)
    a[i] <- coef(modelo)[1]
    b[i] <- coef(modelo)[2]
    se[i] <- summary.lm(modelo)$coefficients[2, 2]
    ms_e[i] <- anova(modelo)[2, 3]
    ms_gxe[i] <- sum((int.mean[i, ] - geno.mean[i] - env.mean + overall.mean)^2) / (env.num - 1)
  }
  stability_geno <- cbind(b, se, ms_e, ms_gxe)
  row.names(stability_geno) <- levels(data[, geno])
  names(a) <- levels(data[, geno])
  names(b) <- levels(data[, geno])
  if (env.num > 2) {
    x <- NULL
    ypred <- NULL
    ymean <- NULL
    for (i in 1:length(data[, trait])) {
      x[i] <- env.mean[names(env.mean) == data[i, env]]
      ypred[i] <- a[names(a) == data[i, geno]] + b[names(b) == data[i, geno]] * x[i]
      ymean[i] <- int.mean[row.names(int.mean) == data[i, geno], colnames(int.mean) == data[i, env]]
    }
    drg_sc <- sum((ypred - ymean)^2)
    hrg_gl <- geno.num - 1
    drg_gl <- (geno.num - 1) * (env.num - 1) - hrg_gl
    drg_cm <- drg_sc / drg_gl
    hrg_sc <- at[4, 2] - drg_sc
    hrg_cm <- hrg_sc / hrg_gl
    hrg_f <- hrg_cm / drg_cm
    hrg_p <- pf(hrg_f, hrg_gl, drg_gl, lower.tail = FALSE)
    drg_f <- drg_cm / at[5, 3]
    drg_p <- pf(drg_f, drg_gl, at[5, 1], lower.tail = FALSE)
  } else {
    drg_sc <- NA
    hrg_gl <- NA
    drg_gl <- NA
    drg_cm <- NA
    hrg_sc <- NA
    hrg_cm <- NA
    hrg_f <- NA
    hrg_p <- NA
    drg_f <- NA
    drg_p <- NA
  }

  # Regression-stability for environments

  a <- NULL        # linear regression intercept
  b <- NULL        # linear regression slope
  se <- NULL       # slope standard error
  ms_e <- NULL     # error mean square
  ms_gxe <- NULL   # gxe variance
  
  for (i in 1:env.num) {
    modelo <- lm(int.mean[, i] ~ geno.mean)
    a[i] <- coef(modelo)[1]
    b[i] <- coef(modelo)[2]
    se[i] <- summary.lm(modelo)$coefficients[2, 2]
    ms_e[i] <- anova(modelo)[2, 3]
    ms_gxe[i] <- sum((int.mean[, i] - env.mean[i] - geno.mean + overall.mean)^2) / (geno.num - 1)
  }
  stability_env <- cbind(b, se, ms_e, ms_gxe)
  row.names(stability_env) <- levels(data[, env])
  names(a) <- levels(data[, env])
  names(b) <- levels(data[, env])
  if (geno.num > 2) {
    x <- NULL
    ypred <- NULL
    ymean <- NULL
    for (i in 1:length(data[, trait])) {
      x[i] <- geno.mean[names(geno.mean) == data[i, geno]]
      ypred[i] <- a[names(a) == data[i, env]] + b[names(b) == data[i, env]] * x[i]
      ymean[i] <- int.mean[row.names(int.mean) == data[i, geno], colnames(int.mean) == data[i, env]]
    }
    dre_sc <- sum((ypred - ymean)^2)
    hre_gl <- env.num - 1
    dre_gl <- (geno.num - 1) * (env.num - 1) - hre_gl
    dre_cm <- dre_sc / dre_gl
    hre_sc <- at[4, 2] - dre_sc
    hre_cm <- hre_sc / hre_gl
    hre_f <- hre_cm / dre_cm
    hre_p <- pf(hre_f, hre_gl, dre_gl, lower.tail = FALSE)
    dre_f <- dre_cm / at[5, 3]
    dre_p <- pf(dre_f, dre_gl, at[5, 1], lower.tail = FALSE)
  } else {
    dre_sc <- NA
    hre_gl <- NA
    dre_gl <- NA
    dre_cm <- NA
    hre_sc <- NA
    hre_cm <- NA
    hre_f <- NA
    hre_p <- NA
    dre_f <- NA
    dre_p <- NA
  }

  # ANOVA plus regression stability

  fileaux <- at[5, ]
  at[5, ] <- c(hrg_gl, hrg_sc, hrg_cm, hrg_f, hrg_p)
  at[6, ] <- c(drg_gl, drg_sc, drg_cm, drg_f, drg_p)
  at[7, ] <- c(hre_gl, hre_sc, hre_cm, hre_f, hre_p)
  at[8, ] <- c(dre_gl, dre_sc, dre_cm, dre_f, dre_p)
  at[9, ] <- fileaux
  row.names(at) <- c("G", "E", "R:E", "GxE", "- Het.Regr.G", "- Dev.Regr.G",
                     "- Het.Regr.E", "- Dev.Regr.E", "Residuals")
  cv <- sqrt(at[5, 3]) / abs(overall.mean) * 100

  # Return

  list(ANOVA = at, CV = cv, Stability_for_genotypes = stability_geno,
       Stability_for_environments = stability_env)
}
