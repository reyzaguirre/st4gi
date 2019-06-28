#' AMMI or GGE biplots
#'
#' This function produces AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000) biplots
#' using the base system and ggplot2.
#' @param x An object of class \code{ammi}.
#' @param bp.type Choose 1 for the trait-PC1 biplot and 2 for the PC1-PC2 biplot.
#' @param bp1.type Choose "effects" or "means" for biplot-1.
#' @param graph.type \code{"base"} or \code{"ggplot"}.
#' @param color Color for lines, symbols and/or labels for environments, genotypes and axes
#' (Only for the base system plot).
#' @param ... Additional plot arguments (Only for the base system plot).
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
#' plot(model.ammi, graph.type = "ggplot")
#' @importFrom graphics abline text
#' @export

plot.st4gi_ammi <- function(x, bp.type = 2, bp1.type = c("effects", "means"),
                            graph.type = c("base", "ggplot"),
                            color = c("darkorange", "black", "gray"), ...) {
  
  # match arguments
  
  bp1.type <- match.arg(bp1.type)
  graph.type <- match.arg(graph.type)
  
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
  geno.num <- x$Number_of_genotypes
  env.num <- x$Number_of_environments
  
  # Local variables (Only for ggplot)
  
  if (graph.type == "ggplot") {
    xx <- NULL
    x0 <- NULL
    xe <- NULL
    yy <- NULL
    y0 <- NULL
    ye <- NULL
    type <- NULL
  } 
  
  #  Biplot-1
  
  if (bp.type == 1) {
    
    main <- paste0(method, " biplot-1 for ", trait.name)
    
    if (bp1.type == "effects") {

      if (graph.type == "base") {
        minx <- min(c(env.mean - overall.mean, geno.mean - overall.mean)) * 1.1
        maxx <- max(c(env.mean - overall.mean, geno.mean - overall.mean)) * 1.1
        limx <- c(minx, maxx)
      }
      
      xlab <- "Genotype and environment effects"
      xcorg <- geno.mean - overall.mean
      xcore <- env.mean - overall.mean
      xline <- 0

    }
    
    if (bp1.type == "means") {
      
      if (graph.type == "base") {
        limx <- range(c(env.mean, geno.mean))
        limx <- limx + c(-max(abs(limx)), max(abs(limx))) * 0.05
      }
      
      xlab <- "Genotype and environment means"
      xcorg <- geno.mean
      xcore <- env.mean
      xline <- overall.mean

    }
    
    if (graph.type == "base") {
      
      limy <- range(c(E[, 1], G[, 1]))
      
      plot(1, type = "n", xlim = limx, ylim = limy, main = main, xlab = xlab,
           ylab = paste("PC1 (", format(PC.cont[1], digits = 3), "%)"))

      points(xcorg, G[, 1], col = color[2], pch = 17)
      text(xcorg, G[, 1], labels = rownames(int.mean), col = color[2], pos = 1, offset = 0.3)
      
      points(xcore, E[, 1], col = color[1], pch = 15)
      text(xcore, E[, 1], labels = colnames(int.mean), col = color[1], pos = 1, offset = 0.3)
      
      abline(h = 0, v = xline, col = color[3], lty = 5)

    } else {
        
      dfr <- data.frame(xx = c(xcorg, xcore), yy = c(G[, 1], E[, 1]),
                        type = c(rep("geno", geno.num), rep("env", env.num)))
      
      gg <- ggplot2::ggplot(dfr) +
        ggplot2::geom_point(ggplot2::aes(x = xx, y = yy, colour = type)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(main) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(paste("PC1 (", format(PC.cont[1], digits = 3), "%)")) +
        ggplot2::geom_vline(xintercept = xline, col = "gray", linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, col = "gray", linetype = 2) +
        ggrepel::geom_text_repel(ggplot2::aes(x = xx, y = yy, label = rownames(dfr), color = type))
      
    }
    
  }
  
  # Biplot-2
  
  if (bp.type == 2) {
    
    main <- paste0(method, " biplot-2 for ", trait.name)
    
    if (graph.type == "base") {

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

    } else {
    
      dfr <- data.frame(xx = c(G[, 1], E[, 1]), yy = c(G[, 2], E[, 2]),
                        type = c(rep("geno", geno.num), rep("env", env.num)))
      
      dfr2 <- data.frame(x0 = rep(0, env.num), y0 = rep(0, env.num), xe = E[, 1], ye = E[, 2])
      
      gg <- ggplot2::ggplot(dfr) +
        ggplot2::geom_point(ggplot2::aes(x = xx, y = yy, colour = type)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(main) +
        ggplot2::xlab(paste("PC1 (", format(PC.cont[1], digits = 3), "%)")) +
        ggplot2::ylab(paste("PC2 (", format(PC.cont[2], digits = 3), "%)")) +
        ggplot2::geom_vline(xintercept = 0, col = "gray", linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, col = "gray", linetype = 2) +
        ggrepel::geom_text_repel(ggplot2::aes(x = xx, y = yy, label = rownames(dfr), color = type)) +
        ggplot2::geom_segment(data = dfr2, ggplot2::aes(x = x0, y = y0, xend = xe, yend = ye), linetype = 2)
      
    }
    
  }
  
  if (graph.type == "ggplot")
    gg
  
}
