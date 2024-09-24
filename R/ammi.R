#' AMMI or GGE with data at plot level
#'
#' This function runs AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000)
#' with data at plot level.
#' @param y The name of the column for the variable to analyze.
#' @param geno The name of the column that identifies the genotypes.
#' @param env The name of the column that identifies the environments.
#' @param rep The name of the column that identifies the replications or blocks. A RCBD is assumed.
#' @param dfr The name of the data frame.
#' @param method \code{"ammi"} or \code{"gge"}.
#' @param f Scaling factor, defaults to 0.5.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @details Significance of PCs are evaluated only with \code{method = "ammi"}
#' and if the data are balanced (with missing values estimated if possible).
#' @return It returns an object of class \code{ammi} with the overall, genotype,
#' environment and interaction means, the interaction effects matrix, the
#' first and second PC values for genotypes and environments, the ANOVA table,
#' and a table with the contribution of each PC. Significance of PCs are
#' included only if \code{method = "ammi"} and the data are balanced
#' (with missing values estimated if possible).
#' @author Raul Eyzaguirre.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based
#' on the GGE biplot, Crop Sci., Vol 40: 597-605.
#' @seealso \code{svd}
#' @examples
#' model.ammi <- ammi("y", "geno", "env", "rep", met8x12)
#' model.ammi
#' model.gge <- ammi("y", "geno", "env", "rep", met8x12, method = "gge")
#' model.gge
#' @importFrom stats aov deviance
#' @export

ammi <- function(y, geno, env, rep, dfr, method = c("ammi", "gge"),
                 f = 0.5, maxp = 0.1) {

  # Match arguments
  
  method <- match.arg(method)

  # Everything as character

  dfr[, geno] <- as.character(dfr[, geno])
  dfr[, env] <- as.character(dfr[, env])
  dfr[, rep] <- as.character(dfr[, rep])

  # Check data

  lc <- ck.f(y, c(geno, env), rep, dfr)

  # Error messages

  if (lc$nt.0 > 0)
    stop("Some GxE cells have zero frequency.")

  if (lc$nrep == 1)
    warning("There is only one replication. Inference is not possible with one replication.")
  
  if (lc$nl[1] < 2 | lc$nl[2] < 2)
    stop("This is not a MET experiment.")

  if (lc$nl[1] < 3 | lc$nl[2] < 3)
    stop("You need at least 3 genotypes and 3 environments to run AMMI or GGE.")

  # Compute ANOVA

  if (lc$nrep > 1 & lc$nt.mult == 0) {
    aov.model <- aov.met(y, geno, env, rep, dfr, maxp)
    if (lc$nmis > 0) {
      y.est <- paste0(y, ".est")
      dfr[, y] <- mve.met(y, geno, env, rep, dfr, maxp)[, y.est]
    }
  } else {
    lc$nrep <- NULL
    aov.model <- NULL
  }

  # Compute interaction means matrix
  
  int.mean <- tapply(dfr[, y], list(dfr[, geno], dfr[, env]), mean, na.rm = TRUE)

  # Run ammi.gxe

  ammi.gxe(int.mean, y, method, f, aov.model, lc$nrep)
  
}
