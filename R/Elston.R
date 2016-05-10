#'  Elston Index
#'
#' Function to compute the Elston index (Elston, R. C., 1963).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param data The name of the data frame containing the data.
#' @param means The genotypic means to compute the index, \code{"single"} or
#' \code{"fitted"}. The default is \code{"single"}. See details for more information.
#' @param model Type of model to fit means, \code{"gxe"} for a model with gxe interaction
#' or \code{"g+e"} for a model without interaction. The default is \code{"gxe"}.
#' See details for more information.
#' @param lb Lower bound. \code{1} for \eqn{k = min(x)} and \code{2} for
#' \eqn{k = (n \times min(x) - max(x)) / (n - 1)}
#' @author Raul Eyzaguirre
#' @details The Elston index is a weight free index. It is assumed that all the
#' traits are in the same direction where the highest the value the better the
#' genotype. To include any trait with an opposite direction it must be transformed
#' by multiplication by \code{-1} before.
#'
#' If \code{means = "fitted"} then a linear model is fitted and used to estimate
#' the means for each genotype. In this case you can choose between a model with
#' interaction (\code{model = "gxe"}) and without interaction (\code{model = "g+e"}).
#' If \code{means = "fitted"} and \code{model = "gxe"} then the arguments \code{env}
#' and \code{rep} must be specified.
#' If \code{means = "fitted"} and \code{model = "g+e"} then the argument \code{env}
#' must be specified.
#' If \code{means = "single"} and \code{env} is specified, then single arithmetic
#' means are computed over the replications for each genotype at each environment
#' and then for each genotype over environments.
#' If \code{means = "single"} and \code{env} is not specified, then single arithmetic
#' means are computed over all the observations for each genotype.
#' @return It returns a data frame with the genotypic means for each trait,
#' the Elston index, and the rank for each genotype according to the index.
#' @references
#' Elston, R. C. (1963). A weight-free index for the purpose of ranking or selection
#' with respect to several traits at a time. Biometrics. 19(1): 85-97.
#' @examples
#' elston(c("rytha", "bc", "dm", "star", "nocr"), "geno", data = spg)
#' @export

elston <- function(traits, geno, env = NULL, rep = NULL, data,
                   means = c("single", "fitted"), model = c("gxe", "g+e"),
                   lb = 1) {

  # match arguments

  means <- match.arg(means)
  model <- match.arg(model)

  # Error messages

  if (means == "fitted" & model == "gxe" & (is.null(env) | is.null(rep)))
    stop("For 'fitted' means with the 'gxe' model you must specify the arguments 'env' and 'rep'.")

  if (means == "fitted" & model == "g+e" & is.null(env))
    stop("For 'fitted' means with the 'g+e' model you must specify the argument 'env'.")

  # inits

  nt <- length(traits) # number of traits
  k <- NULL
  ng <- nlevels(factor(data[, geno])) # number of genotypes

  # compute means

  outind <- data.frame(geno = levels(factor(data[, geno])))
  colnames(outind) <- geno

  if (means == "single" & is.null(env)) {
    outind <- docomp("mean", traits, geno, data = data)
    colnames(outind) <- c("geno", paste("m", traits, sep = "."))
  }

  if (means == "single" & !is.null(env)) {
    outind <- docomp("mean", traits, c(geno, env), data = data)
    outind <- docomp("mean", traits, geno, data = outind)
    colnames(outind) <- c("geno", paste("m", traits, sep = "."))
  }

  if (means == "fitted") {
    for (i in 1:nt) {
      if (model == "gxe") {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", geno, ":", env,
                               ") + (1|", env, "/", rep, ")"))
        fm <- lme4::lmer(ff, data = data)
      }
      if (model == "g+e") {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", env, ")"))
        fm <- lme4::lmer(ff, data = data)
      }
      temp <- as.data.frame(lme4::fixef(fm))
      colnames(temp) <- paste("f", traits[i], sep = ".")
      temp[, geno] <- substring(rownames(temp), nchar(geno) + 1)
      outind <- merge(outind, temp, all = TRUE)
    }
  }

  # Standardized means

  for (i in 2:(1 + nt))
    outind[, i + nt] <- (outind[, i] - mean(outind[, i], na.rm = TRUE)) / sd(outind[, i], na.rm = TRUE)
  colnames(outind)[(2 + nt):(1 + 2 * nt)] <- c(paste("s", traits, sep = "."))

  # compute lower bounds

  if (lb == 1)
    for (i in 1:nt)
      k[i] <- min(outind[, 1 + nt + i], na.rm = TRUE)

  if (lb == 2)
    for (i in 1:nt)
      k[i] <- (ng * min(outind[, 1 + nt + i], na.rm = TRUE) -
                 max(outind[, 1 + nt + i], na.rm = TRUE)) / (ng - 1)

  # Elston index

  outind$E.Index <- outind[, nt + 2] - k[1]
  if (nt > 1)
    for (i in 2:nt)
      outind$E.Index <- outind$E.Index * (outind[, 1 + nt + i] - k[i])

  outind <- outind[, c(1:(1 + nt), 2 + 2 * nt)]
  outind$E.Rank <- rank(-outind$E.Index, na.last = "keep")

  # results

  outind
}

