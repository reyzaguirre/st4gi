#' Elston Index
#'
#' Function to compute the Elston index (Elston, R. C., 1963).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @param means The genotypic means to compute the index, \code{"single"} or
#' \code{"fitted"}. The default is \code{"single"}. See details for more information.
#' @param lb Lower bound. \code{1} for \eqn{k = min(x)} and \code{2} for
#' \eqn{k = (n \times min(x) - max(x)) / (n - 1)}
#' @details The Elston index is a weight free index. It is assumed that all the
#' traits are in the same direction where the highest the value the better the
#' genotype. To include any trait with an opposite direction it must be transformed
#' by multiplication by \code{-1} before.
#'
#' If \code{means = "fitted"} then a linear model is fitted and used to estimate
#' the means for each genotype.
#' If \code{means = "single"} and \code{env} is specified, then single arithmetic
#' means are computed over the replications for each genotype at each environment
#' and then for each genotype over environments.
#' If \code{means = "single"} and \code{env} is not specified, then single arithmetic
#' means are computed over all the observations for each genotype.
#' @return It returns a data frame with the genotypic means for each trait,
#' the Elston index, and the rank for each genotype according to the index.
#' @author Raul Eyzaguirre.
#' @references
#' Elston, R. C. (1963). A weight-free index for the purpose of ranking or selection
#' with respect to several traits at a time. Biometrics. 19(1): 85-97.
#' @examples
#' elston(c("rytha", "bc", "dm", "star", "nocr"), "geno", dfr = spg)
#' @importFrom stats as.formula sd
#' @export

elston <- function(traits, geno, env = NULL, rep = NULL, dfr,
                   means = c("single", "fitted"), lb = 1) {

  # Match arguments

  means <- match.arg(means)
  
  # As character
  
  dfr[, geno] <- as.character(dfr[, geno])

  # inits

  nt <- length(traits) # number of traits
  k <- NULL
  ng <- length(unique(dfr[, geno])) # number of genotypes

  # compute means
  
  if (means == "fitted" & is.null(env) & is.null(rep))
    means <- "single"

  dfr.out <- data.frame(geno = unique(dfr[, geno]))
  colnames(dfr.out) <- geno

  if (means == "single" & is.null(env)) {
    dfr.out <- docomp("mean", traits, geno, dfr = dfr)
    colnames(dfr.out) <- c("geno", paste("m", traits, sep = "."))
  }

  if (means == "single" & !is.null(env)) {
    dfr.out <- docomp("mean", traits, c(geno, env), dfr = dfr)
    dfr.out <- docomp("mean", traits, geno, dfr = dfr.out)
    colnames(dfr.out) <- c("geno", paste("m", traits, sep = "."))
  }

  if (means == "fitted") {
    for (i in 1:nt) {
      if (!is.null(env) & !is.null(rep)) {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", geno, ":", env,
                               ") + (1|", env, "/", rep, ")"))
        fm <- lme4::lmer(ff, dfr)
      }
      if (is.null(env) & !is.null(rep)) {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", rep, ")"))
        fm <- lme4::lmer(ff, dfr)
      }
      temp <- as.data.frame(lme4::fixef(fm))
      colnames(temp) <- paste("f", traits[i], sep = ".")
      temp[, geno] <- substring(rownames(temp), nchar(geno) + 1)
      dfr.out <- merge(dfr.out, temp, all = TRUE)
    }
  }

  # Standardized means

  for (i in 2:(1 + nt))
    dfr.out[, i + nt] <- (dfr.out[, i] - mean(dfr.out[, i], na.rm = TRUE)) / sd(dfr.out[, i], na.rm = TRUE)
  colnames(dfr.out)[(2 + nt):(1 + 2 * nt)] <- c(paste("s", traits, sep = "."))

  # compute lower bounds

  if (lb == 1)
    for (i in 1:nt)
      k[i] <- min(dfr.out[, 1 + nt + i], na.rm = TRUE)

  if (lb == 2)
    for (i in 1:nt)
      k[i] <- (ng * min(dfr.out[, 1 + nt + i], na.rm = TRUE) -
                 max(dfr.out[, 1 + nt + i], na.rm = TRUE)) / (ng - 1)

  # Elston index

  dfr.out$E.Index <- dfr.out[, nt + 2] - k[1]
  if (nt > 1)
    for (i in 2:nt)
      dfr.out$E.Index <- dfr.out$E.Index * (dfr.out[, 1 + nt + i] - k[i])

  dfr.out <- dfr.out[, c(1:(1 + nt), 2 + 2 * nt)]
  dfr.out$E.Rank <- rank(-dfr.out$E.Index, na.last = "keep")

  # results

  dfr.out
  
}
