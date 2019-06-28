#' AMMI or GGE with data from an interaction means matrix
#'
#' This function runs AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000)
#' with data from an interaction means matrix.
#' @param int.mean GxE means matrix, genotypes in rows, environments in columns.
#' @param trait.name Name of the trait.
#' @param nrep Number of replications.
#' @param rdf Residual degrees of freedom.
#' @param rms Residual mean square.
#' @param method \code{"ammi"} or \code{"gge"}.
#' @param f Scaling factor, defaults to 0.5.
#' @details Significance of PCs are evaluated only with \code{method = "ammi"} and if
#' \code{nrep}, \code{rms} and \code{rdf} are specified.
#' @return It returns an object of class ammi with the overall, genotype,
#' environment and interaction means, the interaction effects matrix, the
#' first and second PC values for genotypes and environments, and a table
#' with the contribution of each PC. Significance of PCs are included in the
#' contributions table only if \code{method = "ammi"} and \code{nrep}, \code{rms}
#' and \code{rdf} are specified.
#' @author Raul Eyzaguirre.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based on the GGE
#' biplot, Crop Sci., Vol 40: 597-605.
#' @seealso \code{svd}
#' @examples
#' # Compute GxE means
#' int.mean <- tapply(met8x12$y, list(met8x12$geno, met8x12$env), mean, na.rm = TRUE)
#' # Run AMMI with GxE means matrix
#' model.ammi <- ammi.gxe(int.mean, trait.name = "y")
#' model.ammi
#' # Run GGE with GxE means matrix
#' model.gge <- ammi.gxe(int.mean, trait.name = "y", method = "gge")
#' model.gge
#' @importFrom stats pf
#' @export

ammi.gxe <- function(int.mean, trait.name = NULL, nrep = NULL, rdf = NULL,
                     rms = NULL, method = c("ammi", "gge"), f = 0.5) {

  # Match arguments
  
  method <- match.arg(method)

  # Data

  overall.mean <- mean(int.mean)
  env.mean <- apply(int.mean, 2, mean)
  geno.mean <- apply(int.mean, 1, mean)
  env.num <- length(env.mean)
  geno.num <- length(geno.mean)
  int.eff <- int.mean + overall.mean - geno.mean
  int.eff <- t(t(int.eff) - env.mean)
  
  if (method == "ammi")
    svd.mat <- int.eff
  
  if (method == "gge") {
    svd.mat <- int.mean
    svd.mat <- t(t(svd.mat) - env.mean)
  }

  # SVD

  PC <- min(env.num, geno.num) - 1
  dec <- svd(svd.mat, nu = PC, nv = PC)
  D <- diag(dec$d[1:PC])
  G <- dec$u %*% (D^f)
  E <- dec$v %*% (D^(1 - f))
  PC.geno <- cbind(G[, 1], G[, 2])
  dimnames(PC.geno) <- list(rownames(int.mean), c("PC1", "PC2"))
  PC.env <- cbind(E[, 1], E[, 2])
  dimnames(PC.env) <- list(colnames(int.mean), c("PC1", "PC2"))
  PC.num <- paste0("PC", c(1:PC))
  PC.sv <- dec$d[1:PC]^2

  # Contribution of PCs

  PC.cont <- PC.sv / sum(PC.sv) * 100
  PC.acum <- cumsum(PC.cont)
  tablaPC <- data.frame(PC = PC.num, SV = PC.sv, Cont = PC.cont, CumCont = PC.acum)

  # Significance of PCs, only for AMMI and if nrep, rms and rdf are known

  if (method == "ammi") {
    if (!is.null(nrep)) {
      int.SS <- (t(as.vector(svd.mat)) %*% as.vector(svd.mat)) * nrep
      PC.SS <- (dec$d[1:PC]^2) * nrep
      PC.DF <- env.num + geno.num - 1 - 2 * c(1:PC)
      MS <- PC.SS / PC.DF
    }
    if (!is.null(rms) & !is.null(rdf)) {
      f <- MS / rms
      probab <- pf(f, PC.DF, rdf, lower.tail = FALSE)
      rowlab <- PC.num
      tablaPC <- cbind(tablaPC, PC.DF, PC.SS, MS, f, probab)
      colnames(tablaPC)[5:9] <- c("df", "SumSq", "MeanSq", "Fvalue", "Pr(>F)")
    }
  }

  # Output

  output <- list(Method = method, Trait = trait.name,
                 Number_of_genotypes = geno.num,
                 Number_of_environments = env.num,
                 Overall_mean = overall.mean, Genotype_means = geno.mean,
                 Environment_means = env.mean, Interaction_means = int.mean,
                 Interaction_effects = int.eff, PC_values_genotypes = PC.geno,
                 PC_values_environments = PC.env, Contribution_PCs = tablaPC)
  
  class(output) <- "st4gi_ammi"
  invisible(output)
  
}
