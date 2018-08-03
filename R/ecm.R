#' Estimate Genotypic and Phenotypic Covariance and Correlation Matrices
#'
#' This function estimates the genotypic and phenotypic covariance and correlation
#' matrices with data from a RCBD with one or several environments.
#' @param traits The traits to include.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param dfr The name of the data frame.
#' @param method The method to compute genotypic covariances. See details.
#' @details If \code{env = NULL} a RCBD with one environment is considered.
#' If \code{method = 1} the covariances between each pair of traits are computed
#' using the variances of each trait and the variances of the sum. If \code{method = 2}
#' the covariance matrices are approximated using the average of the correlation
#' matrices computed with each replication in each environment.
#' @return It returns the genotypic and phenotypic covariance and correlation matrices.
#' @author Raul Eyzaguirre.
#' @examples
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' ecm(traits, "geno", "loc", "rep", spg)
#' @export

ecm <- function(traits, geno, env = NULL, rep, dfr, method = 1) {

  # Everything as character

  g <- as.character(dfr[, geno])
  if (!is.null(env))
    e <- as.character(dfr[, env])
  r <- as.character(dfr[, rep])

  # Inits
  
  nt <- length(traits) # number of traits
  G <- matrix(nrow = nt, ncol = nt) # genotypic covariance matrix
  P <- matrix(nrow = nt, ncol = nt) # phenotypic covariance matrix
  ng <- length(unique(g)) # number of genotypes
  if (!is.null(env))
    ne <- length(unique(e)) # number of environments
  nrep <- length(unique(r)) # number of replications in each environment

  # Fitted models by REML for variance components
  
  if (!is.null(env)) {
    for (i in 1:nt) {
      y <- dfr[, traits[i]]
      fm <- lme4::lmer(y ~ (1|g) + (1|g:e) + (1|e/r))
      vc <- lme4::VarCorr(fm)
      G[i, i] <- vc$g[1]
      P[i, i] <- vc$g[1] + vc$e[1] / ne + attr(vc, "sc")^2 / ne / nrep
    }
  }
  if (is.null(env)) {
    for (i in 1:nt) {
      y <- dfr[, traits[i]]
      fm <- lme4::lmer(y ~ (1|g) + (1|r))
      vc <- lme4::VarCorr(fm)
      G[i, i] <- vc$g[1]
      P[i, i] <- vc$g[1] + attr(vc, "sc")^2 / nrep
    }
  }
  
  # Method based on Z = Y1 + Y2

  if (method == 1) {
    if (!is.null(env)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- suma(dfr[, traits[i]], dfr[, traits[j]])
          fm <- lme4::lmer(z ~ (1|g) + (1|g:e) + (1|e/r))
          vcz <- lme4::VarCorr(fm) # variance components for z = x + y
          G[i, j] <- G[j, i] <- (vcz$g[1] - G[i, i] - G[j, j]) / 2
          P[i, j] <- P[j, i] <- (vcz$g[1] + vcz$e[1] / ne + attr(vcz, "sc")^2 / ne / nrep -
                                   P[i, i] - P[j, j]) / 2
        }
      }
    }
    if (is.null(env)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- suma(dfr[, traits[i]], dfr[, traits[j]])
          fm <- lme4::lmer(z ~ (1|g) + (1|r))
          vcz <- lme4::VarCorr(fm) # variance components for z = x + y
          G[i, j] <- G[j, i] <- (vcz$g[1] - G[i, i] - G[j, j]) / 2
          P[i, j] <- P[j, i] <- (vcz$g[1] + attr(vcz, "sc")^2 / nrep - P[i, i] - P[j, j]) / 2
        }
      }
    }
    d1 <- diag(diag(G)^{-0.5}, nt, nt)
    d2 <- diag(diag(P)^{-0.5}, nt, nt)
    GC <- d1 %*% G %*% d1 # Genotypic correlation matrix
    PC <- d2 %*% P %*% d2 # Phenotypic correlation matrix
  }
  
  # Method based on average correlation matrix
  
  if (method == 2) {
    if (!is.null(env)) {
      df <- dfr[, c(sapply(traits, c), env, rep)] # data frames for each env and rep
      df <- split(df, factor(paste(dfr[, env], dfr[, rep]))) # split by env and rep
    }
    if (is.null(env)) {
      df <- dfr[, c(sapply(traits, c), rep)] # data frames for each env and rep
      df <- split(df, dfr[, rep]) # split by rep
    }
    ner <- length(df)
    cl <- list() # correlation list
    for (i in 1:ner)
      cl[[i]] <- cor(df[[i]][, 1:nt], use = "pairwise.complete.obs")
    GC <- PC <- apply(simplify2array(cl), 1:2, mean, na.rm = TRUE) # Correlation matrices
    d1 <- diag(diag(G)^0.5, nt, nt)
    d2 <- diag(diag(P)^0.5, nt, nt)
    G <- d1 %*% GC %*% d1 # Genotypic covariance matrix
    P <- d2 %*% PC %*% d2 # Phenotypic covariance matrix
  }
  
  dimnames(G) <- dimnames(P) <- dimnames(GC) <- dimnames(PC) <- list(traits, traits)

  # results
  
  list(G.Cov = G, P.Cov = P, G.Cor = GC, P.Cor = PC)
  
}