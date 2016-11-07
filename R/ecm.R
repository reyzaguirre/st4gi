#' Estimate Genotypic and Phenotypic Covariance and Correlation Matrices
#'
#' This function estimates the genotypic and phenotypic covariance and correlation
#' matrices with data from a RCBD with one or several environments.
#' @param traits The traits to include.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame containing the data.
#' @param method The method to compute genotypic covariances. See details.
#' @author Raul Eyzaguirre.
#' @details If \code{env = NULL} a RCBD with one environment is considered.
#' If \code{method = 1} the covariance between each pair of traits is computed using the variance
#' of each trait and the variance of the sum. If \code{method = 2} the covariance matrix is
#' approximated using the average of the correlation matrices computed with each replication.
#' @return It returns the genotypic and phenotypic covariance and correlation matrices.
#' @examples
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' ecm(traits, "geno", "loc", "rep", spg)
#' @export

ecm <- function(traits, geno, env = NULL, rep, data, method = 1) {

  # Everything as factor

  g <- factor(data[, geno])
  if (!is.null(env))
    e <- factor(data[, env])
  r <- factor(data[, rep])

  # inits
  
  nt <- length(traits) # number of traits
  G <- matrix(nrow = nt, ncol = nt) # genotypic covariance matrix
  P <- matrix(nrow = nt, ncol = nt) # phenotypic covariance matrix
  ng <- nlevels(g) # number of genotypes
  if (!is.null(env))
    ne <- nlevels(e) # number of environments
  nr <- nlevels(r) # number of replications in each environment

  # fitted models by REML for variance components
  
  if (!is.null(env)) {
    for (i in 1:nt) {
      y <- data[, traits[i]]
      fm <- lme4::lmer(y ~ (1|g) + (1|g:e) + (1|e/r))
      vc <- lme4::VarCorr(fm)
      G[i, i] <- vc$g[1]
      P[i, i] <- vc$g[1] + vc$e[1] / ne + attr(vc, "sc")^2 / ne / nr
    }
  }
  if (is.null(env)) {
    for (i in 1:nt) {
      y <- data[, traits[i]]
      fm <- lme4::lmer(y ~ (1|g) + (1|r))
      vc <- lme4::VarCorr(fm)
      G[i, i] <- vc$g[1]
      P[i, i] <- vc$g[1] + attr(vc, "sc")^2 / nr
    }
  }

  if (method == 1) {
    if (!is.null(env)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- suma(data[, traits[i]], data[, traits[j]])
          fm <- lme4::lmer(z ~ (1|g) + (1|g:e) + (1|e/r))
          vcz <- lme4::VarCorr(fm) # variance components for z = x + y
          G[i, j] <- G[j, i] <- (vcz$g[1] - G[i, i] - G[j, j]) / 2
          P[i, j] <- P[j, i] <- (vcz$g[1] + vcz$e[1] / ne + attr(vcz, "sc")^2 / ne / nr -
                                   P[i, i] - P[j, j]) / 2
        }
      }
    }
    if (is.null(env)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- suma(data[, traits[i]], data[, traits[j]])
          fm <- lme4::lmer(z ~ (1|g) + (1|r))
          vcz <- lme4::VarCorr(fm) # variance components for z = x + y
          G[i, j] <- G[j, i] <- (vcz$g[1] - G[i, i] - G[j, j]) / 2
          P[i, j] <- P[j, i] <- (vcz$g[1] + attr(vcz, "sc")^2 / nr - P[i, i] - P[j, j]) / 2
        }
      }
    }
    d1 <- diag(diag(G)^{-0.5}, nt, nt)
    d2 <- diag(diag(P)^{-0.5}, nt, nt)
    GC <- d1 %*% G %*% d1 # Genotypic correlation matrix
    PC <- d2 %*% P %*% d2 # Phenotypic correlation matrix
  }
  
  # correlation matrix
  
  if (method == 2) {
    if (!is.null(env)) {
      df <- data[, c(sapply(traits, c), env, rep)] # data frames for each env and rep
      df <- split(df, factor(paste(data[, env], data[, rep]))) # split by env and rep
    }
    if (is.null(env)) {
      df <- data[, c(sapply(traits, c), rep)] # data frames for each env and rep
      df <- split(df, data[, rep]) # split by rep
    }
    ner <- length(df)
    cl <- list() # correlation list
    for (i in 1:ner)
      cl[[i]] <- cor(df[[i]][, 1:nt], use = "pairwise.complete.obs")
    GC <- apply(simplify2array(cl), 1:2, mean, na.rm = TRUE) # Genotypic correlation matrix
    S <- diag(diag(G)^0.5, nt, nt)
    G <- S %*% GC %*% S # Genotypic covariance matrix
    pv <- diag(P)
    P <- G # Phenotypic covariance matrix
    diag(P) <- pv
    d2 <- diag(diag(P)^{-0.5}, nt, nt)
    PC <- d2 %*% P %*% d2 # Phenotypic correlation matrix
  }
  
  dimnames(G) <- dimnames(P) <- dimnames(GC) <- dimnames(PC) <- list(traits, traits)

  # results
  
  list(G.Cov = G, P.Cov = P, G.Cor = GC, P.Cor = PC)
}