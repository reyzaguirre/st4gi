#' AMMI or GGE with data at plot level
#'
#' This function runs AMMI (Gollob, H. R., 1968) or GGE biplot (Yan , W. et al., 2000)
#' with data at plot level.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame containing the data.
#' @param method AMMI or GGE.
#' @param f Scaling factor, defaults to 0.5.
#' @param biplot Choose 1 for the trait-PC1 biplot and 2 for the PC1-PC2 biplot.
#' @param biplot1 Choose "effects" or "means" for biplot1.
#' @param title Main title for biplot1 or biplot2.
#' @param xlab Xlab for biplot1.
#' @param color Color for lines, symbols and/or labels for environments, genotypes and axes.
#' @param ... Additional graphic parameters.
#' @author Raul Eyzaguirre
#' @details Significance of PCs are evaluated only with \code{method = "AMMI"} and if
#' the data are balanced.
#' @return It returns the first and second PC values for genotypes and environments,
#' a table with the contribution of each PC, a dispersion plot of means or effects
#' against the first PC, or a dispersion plot of PC1 against PC2.
#' Significance of PCs are included in the contributions table only if
#' method is set to AMMI.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based on the GGE
#' biplot, Crop Sci., Vol 40: 597-605.
#' @seealso \code{svd}
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Run AMMI for trait y, biplot1 by default
#' ammi("y", "geno", "env", "rep", met8x12)
#'
#' # Run AMMI for trait y, biplot2
#' ammi("y", "geno", "env", "rep", met8x12, biplot = 2)
#' @export

ammi <- function(trait, geno, env, rep, data, method = "AMMI", f = .5,
                 biplot = 1, biplot1 = "effects", title = NULL, xlab = NULL,
                 color = c("darkorange", "black", "gray"), ...){

  # Everything as factor

  data[,geno] <- factor(data[,geno])
  data[,env] <- factor(data[,env])
  data[,rep] <- factor(data[,rep])

  # Check data

  lc <- checkdata02(trait, geno, env, data)

  # Error messages

  if (lc$c1 == 0)
    stop("Error: Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")

  if (lc$c1 == 1 & lc$c2 == 0)
    warning("Warning: There is only one replication. Inference is not possible with one replication.")

  if (method == "AMMI" & lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 0)
    warning("Warning: The data set is unbalanced. Significance of PCs is not evaluated.")

  geno.num <- nlevels(data[,geno])
  env.num <- nlevels(data[,env])
  rep.num <- nlevels(data[,rep])

  if (geno.num < 2 | env.num < 2)
    stop(paste("Error: This is not a MET experiment."))

  if (geno.num < 3 | env.num < 3)
    stop(paste("Error: You need at least 3 genotypes and 3 environments to run AMMI or GGE."))

  # Compute interaction means matrix

  int.mean <- tapply(data[,trait], list(data[,geno], data[,env]), mean, na.rm = T)

  # Compute ANOVA

  if (lc$c1 + lc$c2 + lc$c3 == 3){
    model <- aov(data[,trait] ~ data[,geno] + data[,env] +
                   data[,rep] %in% data[,env] + data[,geno]:data[,env])
    rdf <- model$df.residual
    rms <- deviance(model)/rdf
  } else {
    rep.num <- NULL
    rdf <- NULL
    rms <- NULL
  }

  # Run ammigxe

  ammigxe(int.mean, trait = trait, rep.num = rep.num, rdf = rdf, rms = rms,
          method = method, f = f, biplot = biplot, biplot1 = biplot1,
          title = title, xlab = xlab, color = color, ...)
}

#' AMMI or GGE with data from an interaction means matrix
#'
#' This function runs AMMI (Gollob, H. R., 1968) or GGE biplot (Yan , W. et al., 2000)
#' with data from an interaction means matrix.
#' @param int.mean GxE means matrix, genotypes in rows, environments in columns.
#' @param trait Name of the trait.
#' @param rep.num Number of replications.
#' @param rdf Residual degrees of freedom.
#' @param rms Residual mean square.
#' @param method AMMI or GGE.
#' @param f Scaling factor, defaults to 0.5.
#' @param biplot 1 for the trait-PC1 biplot and 2 for the PC1-PC2 biplot.
#' @param biplot1 Choose "effects" or "means" for biplot1.
#' @param title Main title for biplot1 or biplot2.
#' @param xlab Xlab for biplot1.
#' @param color Color for lines, symbols and/or labels for environments, genotypes and axes.
#' @param ... Additional graphic parameters.
#' @author Raul Eyzaguirre
#' @details Significance of PCs are evaluated only with \code{method = "AMMI"} and if
#' \code{rep.num}, \code{rms} and \code{rdf} are specified.
#' @return It returns the first and second PC values for genotypes and environments,
#' a table with the contribution of each PC, a dispersion plot of means or effects
#' against the first PC, or a dispersion plot of PC1 against PC2.
#' Significance of PCs are included in the contributions table only if \code{rep.num},
#' \code{rms} and \code{rdf} are specified and method is set to AMMI.
#' @references
#' Gollob, H. R. (1968). A Statistical Model which combines Features of Factor Analytic
#' and Analysis of Variance Techniques, Psychometrika, Vol 33(1): 73-114.
#'
#' Yan, W. et al. (2000). Cultivar evaluation and mega-environment investigation based on the GGE
#' biplot, Crop Sci., Vol 40: 597-605.
#' @seealso \code{svd}
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Compute GxE means
#' int.mean <- tapply(met8x12$y, list(met8x12$geno, met8x12$env), mean, na.rm = TRUE)
#'
#' # Run AMMI with GxE means matrix, biplot2
#' ammigxe(int.mean, trait = "y", biplot = 2)
#'
#' # Run GGE with GxE means matrix, biplot2
#' ammigxe(int.mean, trait = "y", method = "GGE", biplot = 2)
#' @export

ammigxe <- function(int.mean, trait = NULL, rep.num = NULL, rdf = NULL, rms = NULL,
                    method = "AMMI", f = .5, biplot = 1, biplot1 = "effects",
                    title = NULL, xlab = NULL, color = c("darkorange", "black", "gray"), ...){

  # Data

  overall.mean <- mean(int.mean)
  env.mean <- apply(int.mean, 2, mean)
  geno.mean <- apply(int.mean, 1, mean)
  env.num <- length(env.mean)
  geno.num <- length(geno.mean)

  if (method == "AMMI"){
    svd.mat <- int.mean + overall.mean
    for (i in 1:env.num) svd.mat[,i] <- svd.mat[,i] - geno.mean
    for (i in 1:geno.num) svd.mat[i,] <- svd.mat[i,] - env.mean
  }

  if (method == "GGE"){
    svd.mat <- int.mean
    for (i in 1:geno.num) svd.mat[i,] <- svd.mat[i,] - env.mean
  }

  # SVD

  PC <- min(env.num, geno.num)-1
  dec <- svd(svd.mat, nu = PC, nv = PC)
  D <- diag(dec$d[1:PC])
  G <- dec$u %*% (D^f)
  E <- dec$v %*% (D^(1-f))
  PC.geno <- cbind(G[,1], G[,2])
  dimnames(PC.geno) <- list(rownames(int.mean), c("PC1", "PC2"))
  PC.env <- cbind(E[,1], E[,2])
  dimnames(PC.env) <- list(colnames(int.mean), c("PC1", "PC2"))
  PC.num <- paste("PC", c(1:PC), sep = "")
  PC.sv <- dec$d[1:PC]^2

  # Contribution of PCs

  PC.cont <- PC.sv/sum(PC.sv)*100
  PC.acum <- cumsum(PC.cont)
  tablaPC <- data.frame(PC = PC.num, SV = PC.sv, Cont = PC.cont, CumCont = PC.acum)

  # Significance of PCs, only for AMMI and if rep.num, rms and rdf are known

  if (method == "AMMI"){
    if (is.null(rep.num) == 0){
      int.SS <- (t(as.vector(svd.mat))%*%as.vector(svd.mat))*rep.num
      PC.SS <- (dec$d[1:PC]^2)*rep.num
      PC.DF <- env.num + geno.num - 1 - 2*c(1:PC)
      MS <- PC.SS/PC.DF
    }
    if (is.null(rms) == 0 & is.null(rdf) == 0){
      F <- MS/rms
      probab <- pf(F, PC.DF, rdf, lower.tail = FALSE)
      rowlab <- PC.num
      tablaPC <- cbind(tablaPC, PC.DF, PC.SS, MS, F, probab)
      colnames(tablaPC)[5:9] <- c("df", "SumSq", "MeanSq", "Fvalue", "Pr(>F)")
    }
  }

  #  Biplot 1

  if (biplot == 1){

    if (is.null(title) == 1)
      title = paste(method, " biplot1 for ", trait, sep = "")

    if (biplot1 == "effects"){
      maxx <- max(abs(c(env.mean - overall.mean, geno.mean - overall.mean)))*1.05
      limx <- c(-maxx, maxx)
      if (is.null(xlab) == 1)
        xlab = "Genotype and environment effects"
      xcorg = geno.mean - overall.mean
      xcore = env.mean - overall.mean
      xline = 0
    }
    if (biplot1 == "means"){
      limx <- range(c(env.mean, geno.mean))
      limx <- limx + c(-max(abs(limx)), max(abs(limx)))*.05
      if (is.null(xlab) == 1)
        xlab = "Genotype and environment means"
      xcorg = geno.mean
      xcore = env.mean
      xline = overall.mean
    }

    limy <- c(-max(abs(c(E[,1], G[,1]))), max(abs(c(E[,1], G[,1]))))

    plot(1, type = "n", xlim = limx, ylim = limy, main = title, xlab = xlab,
         ylab = paste("PC1 (",format(PC.cont[1],digits = 3),"%)"), ...)
    points(xcorg, G[,1], col = color[2], pch = 17, ...)
    text(xcorg, G[,1], labels = rownames(int.mean), col = color[2], pos = 1, offset = 0.3)
    points(xcore, E[,1], col = color[1], pch = 15, ...)
    text(xcore, E[,1], labels = colnames(int.mean), col = color[1], pos = 1, offset = .3)
    abline(h = 0, v = xline, col = color[3], lty = 2)
  }

  # Biplot 2

  if (biplot == 2){

    if (is.null(title) == 1)
      title = paste(method, " biplot2 for ", trait, sep = "")

    limx <- range(c(E[,1], G[,1]))
    limx <- limx + c(-max(abs(limx)), max(abs(limx)))*.05
    limy <- range(c(E[,2], G[,2]))

    plot(1, type = "n", xlim = limx, ylim = limy, main = title,
         xlab = paste("PC1 (", format(PC.cont[1], digits = 3), "%)"),
         ylab = paste("PC2 (", format(PC.cont[2], digits = 3), "%)"),
         asp = 1, ...)
    points(G[,1], G[,2], col = color[2], pch = 17, ...)
    text(G[,1], G[,2], labels = rownames(int.mean), col = color[2], pos = 1, offset = .3)
    points(E[,1], E[,2], col = color[1], pch = 15, ...)
    text(E[,1], E[,2], labels = colnames(int.mean), col = color[1], pos = 1, offset = .3)
    abline(h = 0, v = 0, col = color[3], lty = 2)
    for (i in 1:env.num) lines(c(0,E[i,1]), c(0,E[i,2]), col = color[1], lty = 3)
  }

  # Output

  list(PC_values_genotypes = PC.geno, PC_values_environments = PC.env,
       Contribution_PCs = tablaPC)
}
