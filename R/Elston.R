#'  Elston Index
#'
#' Function to compute the Elston index (Elston, R. C., 1963).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks.
#' @param data The name of the data frame containing the data.
#' @param means The genotypic means to compute the index, \code{"single"}
#' or \code{"fitted"}. The default is \code{"single"}. See details for more information.
#' @param model Type of model to fit means if \code{means = "fitted"}, \code{"gxe"} for
#' a model with gxe interaction or \code{"g+e"} for a model without interaction.
#' The default is \code{"gxe"}. See details for more information.
#' @param lb Lower bound. \code{1} for k = min(x) and \code{2} for k = (n * min(x) - max(x)) / (n - 1)
#' @author Raul Eyzaguirre
#' @details The Elston index is a weight free index.
#' 
#' If \code{means = "fitted"} and \code{model = "gxe"} then the arguments \code{env} and
#' \code{rep} must be specified.
#' If \code{means = "fitted"} and \code{model = "g+e"} then only the argument \code{env}
#' must be specified.
#' If \code{means = "single"} and \code{env} and \code{rep} are specified, then
#' single arithmetic means are computed over the replications for each genotype
#' at each location and then for each genotype over locations. In any other case
#' single arithmetic means are computed over all the observations for each genotype.
#' @return It returns a data frame with the genotypic means for each trait, the Elston index,
#' and the rank for each genotype according to the index.
#' @references
#' Elston, R. C. (1963). A weight-free index for the purpose of ranking or selection
#' with respect to several traits at a time. Biometrics. 19(1): 85-97.
#' @examples
#' # The data
#' head(spg)
#' str(spg)
#'
#' # Run Elston index with all the traits
#' elston(c("rytha", "bc", "dm", "star", "nocr"), "geno", data = spg)
#' @export

elston <- function(traits, geno, env = NULL, rep = NULL, data,
                   means = "single", model = "gxe", lb = 1) {

  # inits

  nt <- length(traits) # number of traits
  k <- NULL
  ng <- nlevels(factor(data[, geno])) # number of genotypes

  # compute means
  
  outind <- data.frame(geno = levels(factor(data[, geno])))
  
  if (means == "single" & (is.null(env) | is.null(rep))) {
    m <- matrix(NA, ng, nt)
    for (i in 1:nt)
      m[, i] <- tapply(data[, traits[i]], data[, geno], mean, na.rm = TRUE)
    outind <- cbind(outind, m)
    colnames(outind) <- c("geno", paste("m", traits, sep = "."))
  }
  
  if (means == "single" & !is.null(env) & !is.null(rep)) {
    temp <- domeans(traits, c(geno, env), data = data)
    temp <- domeans(traits, geno, data = temp)
    outind <- merge(outind, temp, all = TRUE)
    colnames(outind) <- c("geno", paste("m", traits, sep = "."))
  }
  
  if (means == "fitted") {
    for (i in 1:nt) {
      abc <- data.frame(c1 = data[, traits[i]], c2 = data[, geno], c3 = data[, env], c4 = data[, rep])
      if (model == "gxe")
        fm <- lme4::lmer(c1 ~ c2 - 1 + (1|c2:c3) + (1|c3 / c4), data = abc)
      if (model == "g+e")
        fm <- lme4::lmer(c1 ~ c2 - 1 + (1|c3), data = abc)
      temp <- as.data.frame(lme4::fixef(fm))
      colnames(temp) <- paste("f", traits[i], sep = ".")
      temp$geno <- substring(rownames(temp), 3)
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

  return(outind)
}
