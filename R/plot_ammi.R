#' AMMI or GGE biplots
#'
#' This function produces AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000) biplots.
#' @param x An object of class \code{ammi}.
#' @param bp.type Choose 1 for the trait-PC1 biplot and 2 for the PC1-PC2 biplot.
#' @param bp1.type Choose "effects" or "means" for biplot-1.
#' @param color Color for lines, symbols and/or labels for environments, genotypes and axes.
#' @param ... Additional plot arguments.
#' @details It produces a biplot for an object of class \code{ammi}. See \code{?ammi}
#' for additional details.
#' @return It returns a dispersion plot of means or effects against the first PC,
#' or a dispersion plot of PC1 against PC2.
#' @author Raul Eyzaguirre.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based
#' on the GGE biplot, Crop Sci., Vol 40: 597-605.
#' @examples
#' model.ammi <- ammi("y", "geno", "env", "rep", met8x12)
#' plot(model.ammi)
#' plot(model.ammi, bp.type = 1)
#' @importFrom graphics abline text
#' @export

plot.ammi <- function(x, bp.type = 2, bp1.type = c("effects", "means"),
                      color = c("darkorange", "black", "gray"), ...) {
  
  # match arguments
  
  bp1.type <- match.arg(bp1.type)
  
  # arguments
  
  method <- x$Method
  trait.name <- x$Trait
  overall.mean <- x$Overall_mean
  geno.mean <- x$Genotype_means
  env.mean <- x$Environment_means
  int.mean <- x$Interaction_means
  G <- x$PC_values_genotypes
  E <- x$PC_values_environments
  PC.cont <- x$Contribution_PCs$Cont
  env.num <- x$Number_of_environments
  
  #  Biplot-1
  
  if (bp.type == 1) {
    
    main <- paste0(method, " biplot-1 for ", trait.name)
    
    if (bp1.type == "effects") {
      minx <- min(c(env.mean - overall.mean, geno.mean - overall.mean)) * 1.1
      maxx <- max(c(env.mean - overall.mean, geno.mean - overall.mean)) * 1.1
      limx <- c(minx, maxx)
      xlab <- "Genotype and environment effects"
      xcorg <- geno.mean - overall.mean
      xcore <- env.mean - overall.mean
      xline <- 0
    }
    if (bp1.type == "means") {
      limx <- range(c(env.mean, geno.mean))
      limx <- limx + c(-max(abs(limx)), max(abs(limx))) * 0.05
      xlab <- "Genotype and environment means"
      xcorg <- geno.mean
      xcore <- env.mean
      xline <- overall.mean
    }
    
    limy <- range(c(E[, 1], G[, 1]))
    
    plot(1, type = "n", xlim = limx, ylim = limy, main = main, xlab = xlab,
         ylab = paste("PC1 (", format(PC.cont[1], digits = 3), "%)"))
    points(xcorg, G[, 1], col = color[2], pch = 17)
    text(xcorg, G[, 1], labels = rownames(int.mean), col = color[2], pos = 1, offset = 0.3)
    points(xcore, E[, 1], col = color[1], pch = 15)
    text(xcore, E[, 1], labels = colnames(int.mean), col = color[1], pos = 1, offset = 0.3)
    abline(h = 0, v = xline, col = color[3], lty = 5)
  }
  
  # Biplot-2
  
  if (bp.type == 2) {
    
    main <- paste0(method, " biplot-2 for ", trait.name)
    
    limx <- range(c(E[, 1], G[, 1]))
    limx <- limx + c(-max(abs(limx)), max(abs(limx))) * 0.1
    limy <- range(c(E[, 2], G[, 2]))
    
    plot(1, type = "n", xlim = limx, ylim = limy, main = main,
         xlab = paste("PC1 (", format(PC.cont[1], digits = 3), "%)"),
         ylab = paste("PC2 (", format(PC.cont[2], digits = 3), "%)"),
         asp = 1)
    points(G[, 1], G[, 2], col = color[2], pch = 17)
    text(G[, 1], G[, 2], labels = rownames(int.mean), col = color[2], pos = 1, offset = 0.3)
    points(E[, 1], E[, 2], col = color[1], pch = 15)
    text(E[, 1], E[, 2], labels = colnames(int.mean), col = color[1], pos = 1, offset = 0.3)
    abline(h = 0, v = 0, col = color[3], lty = 5)
    for (i in 1:env.num)
      lines(c(0, E[i, 1]), c(0, E[i, 2]), col = color[1], lty = 2)
  }
  
}
