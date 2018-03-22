#' AMMI or GGE with data at plot level
#'
#' This function runs AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000)
#' with data at plot level.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame containing the data.
#' @param method \code{"ammi"} or \code{"gge"}.
#' @param f Scaling factor, defaults to 0.5.
#' @author Raul Eyzaguirre.
#' @details Significance of PCs are evaluated only with \code{method = "ammi"} and if
#' the data are balanced.
#' @return It returns an object of class \code{ammi} with the overall, genotype,
#' environment and interaction means, the interaction effects matrix, the
#' first and second PC values for genotypes and environments, and a table
#' with the contribution of each PC. Significance of PCs are included in the
#' contributions table only if \code{method = "ammi"} and the data are balanced.
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

ammi <- function(trait, geno, env, rep, data, method = c("ammi", "gge"), f = 0.5) {

  # match arguments
  
  method <- match.arg(method)

  # Everything as factor

  data[, geno] <- factor(data[, geno])
  data[, env] <- factor(data[, env])
  data[, rep] <- factor(data[, rep])

  # Check data

  lc <- check.2f(trait, geno, env, rep, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")

  if (lc$c2 == 0)
    warning("There is only one replication. Inference is not possible with one replication.")
  
  if (method == "ammi" & (lc$c3 == 0 | lc$c4 == 0))
    warning("The data set is unbalanced. Significance of PCs is not evaluated.")

  if (lc$na < 2 | lc$nb < 2)
    stop(paste("This is not a MET experiment."))

  if (lc$na < 3 | lc$nb < 3)
    stop(paste("You need at least 3 genotypes and 3 environments to run AMMI or GGE."))

  # Compute interaction means matrix

  int.mean <- tapply(data[, trait], list(data[, geno], data[, env]), mean, na.rm = TRUE)

  # Compute ANOVA

  if (lc$c2 + lc$c3 + lc$c4 == 3) {
    model <- aov(data[, trait] ~ data[, geno] + data[, env] +
                   data[, rep] %in% data[, env] + data[, geno]:data[, env])
    rdf <- model$df.residual
    rms <- deviance(model) / rdf
  } else {
    lc$nr <- NULL
    rdf <- NULL
    rms <- NULL
  }

  # Run ammi.gxe

  ammi.gxe(int.mean, trait = trait, nr = lc$nr, rdf = rdf, rms = rms,
           method = method, f = f)
}

#' AMMI or GGE with data from an interaction means matrix
#'
#' This function runs AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000)
#' with data from an interaction means matrix.
#' @param int.mean GxE means matrix, genotypes in rows, environments in columns.
#' @param trait Name of the trait.
#' @param nr Number of replications.
#' @param rdf Residual degrees of freedom.
#' @param rms Residual mean square.
#' @param method \code{"ammi"} or \code{"gge"}.
#' @param f Scaling factor, defaults to 0.5.
#' @author Raul Eyzaguirre.
#' @details Significance of PCs are evaluated only with \code{method = "ammi"} and if
#' \code{nr}, \code{rms} and \code{rdf} are specified.
#' @return It returns an object of class ammi with the overall, genotype,
#' environment and interaction means, the interaction effects matrix, the
#' first and second PC values for genotypes and environments, and a table
#' with the contribution of each PC. Significance of PCs are included in the
#' contributions table only if \code{method = "ammi"} and \code{nr}, \code{rms}
#' and \code{rdf} are specified.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based on the GGE
#' biplot, Crop Sci., Vol 40: 597-605.
#' @seealso \code{svd}
#' @examples
#' ## Compute GxE means
#' int.mean <- tapply(met8x12$y, list(met8x12$geno, met8x12$env), mean, na.rm = TRUE)
#' 
#' ## Run AMMI with GxE means matrix
#' model.ammi <- ammi.gxe(int.mean, trait = "y")
#' model.ammi
#' 
#' ## Run GGE with GxE means matrix
#' model.gge <- ammi.gxe(int.mean, trait = "y", method = "gge")
#' model.gge
#' @importFrom stats pf
#' @export

ammi.gxe <- function(int.mean, trait = NULL, nr = NULL, rdf = NULL, rms = NULL,
                     method = c("ammi", "gge"), f = 0.5) {

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
  PC.num <- paste("PC", c(1:PC), sep = "")
  PC.sv <- dec$d[1:PC]^2

  # Contribution of PCs

  PC.cont <- PC.sv / sum(PC.sv) * 100
  PC.acum <- cumsum(PC.cont)
  tablaPC <- data.frame(PC = PC.num, SV = PC.sv, Cont = PC.cont, CumCont = PC.acum)

  # Significance of PCs, only for AMMI and if nr, rms and rdf are known

  if (method == "ammi") {
    if (!is.null(nr)) {
      int.SS <- (t(as.vector(svd.mat)) %*% as.vector(svd.mat)) * nr
      PC.SS <- (dec$d[1:PC]^2) * nr
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

  output <- list(Method = method, Trait = trait, Number_of_environments = env.num,
                 Overall_mean = overall.mean, Genotype_means = geno.mean,
                 Environment_means = env.mean, Interaction_means = int.mean,
                 Interaction_effects = int.eff, PC_values_genotypes = PC.geno,
                 PC_values_environments = PC.env, Contribution_PCs = tablaPC)
  
  class(output) <- "ammi"
  invisible(output)
}

#' AMMI or GGE biplots
#'
#' This function produces AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000) biplots.
#' @param x An object of class \code{ammi}.
#' @param biplot Choose 1 for the trait-PC1 biplot and 2 for the PC1-PC2 biplot.
#' @param biplot1 Choose "effects" or "means" for biplot1.
#' @param title Main title for biplot1 or biplot2.
#' @param xlab Xlab for biplot1.
#' @param color Color for lines, symbols and/or labels for environments, genotypes and axes.
#' @param size Relative size for symbols and labels.
#' @param ... Additional plot arguments.
#' @author Raul Eyzaguirre.
#' @details It produces a biplot for an object of class \code{ammi}. See \code{?ammi}
#' for additional details.
#' @return It returns a dispersion plot of means or effects against the first PC,
#' or a dispersion plot of PC1 against PC2.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based
#' on the GGE biplot, Crop Sci., Vol 40: 597-605.
#' @examples
#' model.ammi <- ammi("y", "geno", "env", "rep", met8x12)
#' plot(model.ammi)
#' plot(model.ammi, biplot = 1)
#' @importFrom graphics abline text
#' @export

plot.ammi <- function(x, biplot = 2, biplot1 = c("effects", "means"),
                      title = NULL, xlab = NULL, color = c("darkorange", "black", "gray"),
                      size = c(1, 1), ...) {
  
  # match arguments
  
  biplot1 <- match.arg(biplot1)
  
  # arguments
  
  method <- x$Method
  trait <- x$Trait
  overall.mean <- x$Overall_mean
  geno.mean <- x$Genotype_means
  env.mean <- x$Environment_means
  int.mean <- x$Interaction_means
  G <- x$PC_values_genotypes
  E <- x$PC_values_environments
  PC.cont <- x$Contribution_PCs$Cont
  env.num <- x$Number_of_environments
  
  #  Biplot 1
  
  if (biplot == 1) {
    
    if (is.null(title))
      title <- paste(method, " biplot1 for ", trait, sep = "")
    
    if (biplot1 == "effects") {
      minx <- min(c(env.mean - overall.mean, geno.mean - overall.mean)) * 1.1
      maxx <- max(c(env.mean - overall.mean, geno.mean - overall.mean)) * 1.1
      limx <- c(minx, maxx)
      if (is.null(xlab))
        xlab <- "Genotype and environment effects"
      xcorg <- geno.mean - overall.mean
      xcore <- env.mean - overall.mean
      xline <- 0
    }
    if (biplot1 == "means") {
      limx <- range(c(env.mean, geno.mean))
      limx <- limx + c(-max(abs(limx)), max(abs(limx))) * 0.05
      if (is.null(xlab))
        xlab <- "Genotype and environment means"
      xcorg <- geno.mean
      xcore <- env.mean
      xline <- overall.mean
    }
    
    limy <- range(c(E[, 1], G[, 1]))
    
    plot(1, type = "n", xlim = limx, ylim = limy, main = title, xlab = xlab,
         ylab = paste("PC1 (", format(PC.cont[1], digits = 3), "%)"))
    points(xcorg, G[, 1], col = color[2], pch = 17, cex = size[1])
    text(xcorg, G[, 1], labels = rownames(int.mean), col = color[2], pos = 1,
         offset = 0.3, cex = size[2])
    points(xcore, E[, 1], col = color[1], pch = 15, cex = size[1])
    text(xcore, E[, 1], labels = colnames(int.mean), col = color[1], pos = 1,
         offset = 0.3, cex = size[2])
    abline(h = 0, v = xline, col = color[3], lty = 5)
  }
  
  # Biplot 2
  
  if (biplot == 2) {
    
    if (is.null(title))
      title <- paste(method, " biplot2 for ", trait, sep = "")
    
    limx <- range(c(E[, 1], G[, 1]))
    limx <- limx + c(-max(abs(limx)), max(abs(limx))) * 0.1
    limy <- range(c(E[, 2], G[, 2]))
    
    plot(1, type = "n", xlim = limx, ylim = limy, main = title,
         xlab = paste("PC1 (", format(PC.cont[1], digits = 3), "%)"),
         ylab = paste("PC2 (", format(PC.cont[2], digits = 3), "%)"),
         asp = 1)
    points(G[, 1], G[, 2], col = color[2], pch = 17, cex = size[1])
    text(G[, 1], G[, 2], labels = rownames(int.mean), col = color[2], pos = 1,
         offset = 0.3, cex = size[2])
    points(E[, 1], E[, 2], col = color[1], pch = 15, cex = size[1])
    text(E[, 1], E[, 2], labels = colnames(int.mean), col = color[1], pos = 1,
         offset = 0.3, cex = size[2])
    abline(h = 0, v = 0, col = color[3], lty = 5)
    for (i in 1:env.num) lines(c(0, E[i, 1]), c(0, E[i, 2]), col = color[1], lty = 2)
  }
}
