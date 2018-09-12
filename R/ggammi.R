#' AMMI ggplot
#' 
#' This function produces AMMI (Gollob, H. R., 1968) or GGE (Yan , W. et al., 2000) biplots
#' using the ggplot2 package.
#' @param x An object of class \code{ammi}.
#' @param bp.type Choose 1 for the trait-PC1 biplot and 2 for the PC1-PC2 biplot.
#' @param bp1.type Choose "effects" or "means" for biplot-1.
#' @details It produces a ggplot biplot for an object of class \code{ammi}.
#' See \code{?ammi} for additional details.
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
#' ggammi(model.ammi)
#' @export

ggammi <- function(x, bp.type = 2, bp1.type = c("effects", "means")) {
  
  # Match arguments
  
  bp1.type <- match.arg(bp1.type)
  
  # Arguments
  
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
  
  # Local variables
  
  xx <- NULL
  x0 <- NULL
  xe <- NULL
  yy <- NULL
  y0 <- NULL
  ye <- NULL
  type <- NULL
  
  # Biplot-1
  
  if (bp.type == 1) {
    
    main <- paste0(method, " biplot-1 for ", trait.name)
    
    if (bp1.type == "effects") {
      xlab <- "Genotype and environment effects"
      xcorg <- geno.mean - overall.mean
      xcore <- env.mean - overall.mean
      xline <- 0
    }
    if (bp1.type == "means") {
      xlab <- "Genotype and environment means"
      xcorg <- geno.mean
      xcore <- env.mean
      xline <- overall.mean
    }
    
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
  
  # Biplot-2
  
  if (bp.type == 2) {
    
    main <- paste0(method, " biplot-2 for ", trait.name)
    
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
  
  gg
  
}
