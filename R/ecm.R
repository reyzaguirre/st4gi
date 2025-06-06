#' Estimate Genotypic and Phenotypic Covariance and Correlation Matrices
#'
#' This function estimates the genotypic and phenotypic covariance and correlation
#' matrices with data from one or several environments.
#' @param dfr The name of the data frame.
#' @param vars The names of the columns for variables to include.
#' @param geno The name of the column that identifies the genotypes.
#' @param env The name of the column that identifies the environments.
#' @param rep The name of the column that identifies the replications or blocks.
#' @param method The method to compute genotypic covariances and phenotypic
#' variances. See details.
#' @details If \code{env = NULL} data from only one environment is considered.
#' If \code{rep = NULL} data from only one replication is considered without
#' blocking structure.
#' If \code{method = 1} phenotypic covariances are computed from BLUEs and
#' genotypic covariances are computed from BLUPs.
#' If \code{method = 2} the covariances between each pair of variables are computed
#' using the variances of each variable and the variances of the sum.
#' If \code{method = 3} the genotypic covariances are approximated using the average
#' of the correlation matrices computed with each replication in each environment,
#' and the phenotypic covariances are computed pooling all the observed data.
#' @return It returns the genotypic and phenotypic covariance and correlation matrices.
#' @author Raul Eyzaguirre.
#' @examples
#' vars <- c("rytha", "bc", "dm", "star", "nocr")
#' ecm(spg, vars, "geno", "loc", "rep")
#' @importFrom stats cor cov
#' @export

ecm <- function(dfr, vars, geno, env = NULL, rep = NULL, method = 1) {

  # Check there is at least 2 env or reps
  
  if(is.null(env) & is.null(rep))
    stop("There must be at least 2 environments or 2 replications.")
  
  # Everything as character

  g <- as.character(dfr[, geno])
  if (!is.null(env))
    e <- as.character(dfr[, env])
  if (!is.null(rep))
    r <- as.character(dfr[, rep])

  # Inits
  
  nt <- length(vars) # number of variables
  
  G.cov <- P.cov <- matrix(nrow = nt, ncol = nt) # covariance matrices

  ng <- length(unique(g)) # number of genotypes
  
  if (!is.null(env))
    ne <- length(unique(e)) # number of environments
  if (!is.null(rep)) {
    # Check all environments have the same number of replications
    if (!is.null(env) & ne > 1) {
      tmp <- tapply(dfr[, rep], dfr[, env], function(x) length(unique(x)))
      tmp <- length(unique(tmp))
      if (tmp > 1) {
        warning("Number of replications is different among environments; changing to method = 1.")
        method = 1
      }
    }
    if (method != 1)
      nr <- length(unique(r)) # number of replications in each environment
  }
  
  #----------------------------------------------------------------------------
  # Fitted models by REML to get:
  # - BLUEs
  # - BLUPs
  # - Variance components
  #----------------------------------------------------------------------------
  
  dfr.out <- data.frame(geno = unique(dfr[, geno]))
  colnames(dfr.out) <- geno
  vc <- list()
  
  # Models with env and rep

  if (!is.null(env) & !is.null(rep)) {
    fef <- as.formula("y ~ -1 + g + (1|g:e) + (1|e/r)")
    ref <- as.formula("y ~ (1|g) + (1|g:e) + (1|e/r)")
  }
  
  # Models only with env

  if (is.null(rep)) {
    fef <- as.formula("y ~ -1 + g + (1|e)")
    ref <- as.formula("y ~ (1|g) + (1|e)")
  }

  # Models only with rep
  
  if (is.null(env)) {
    fef <- as.formula("y ~ -1 + g + (1|r)")
    ref <- as.formula("y ~ (1|g) + (1|r)")
  }
  
  # Fit models and save blues, blups, and variance components
  
  for (i in 1:nt) {
    
    y <- dfr[, vars[i]]
    
    # Fixed effects model
    
    model <- lme4::lmer(fef)
    tmp <- as.data.frame(lme4::fixef(model))
    colnames(tmp) <- paste("blue", vars[i], sep = ".")
    tmp[, geno] <- substring(rownames(tmp), 2)
    dfr.out <- merge(dfr.out, tmp, all = TRUE)
    
    # Random effects model
    
    model <- lme4::lmer(ref)
    tmp <- lme4::fixef(model) + lme4::ranef(model)$g
    colnames(tmp) <- paste("blup", vars[i], sep = ".")
    tmp[, geno] <- rownames(tmp)
    dfr.out <- merge(dfr.out, tmp, all = TRUE)
    vc[[i]] <- lme4::VarCorr(model)
  }
  
  # Get G.cov and P.cov diagonals for methods 2 and 3
  
  if (method %in% c(2, 3)) {
    
    for (i in 1:nt) {
     
      # Models with env and rep
      
      if (!is.null(env) & !is.null(rep)) {
        G.cov[i, i] <- vc[[i]]$g[1]
        P.cov[i, i] <- vc[[i]]$g[1] + vc[[i]]$e[1] / ne + attr(vc[[i]], "sc")^2 / ne / nr
      }
      
      # Models only with env
      
      if (is.null(rep)) {
        G.cov[i, i] <- vc[[i]]$g[1]
        P.cov[i, i] <- vc[[i]]$g[1] + attr(vc[[i]], "sc")^2 / ne
      }
      
      # Models only with rep
      
      if (is.null(env)) {
        G.cov[i, i] <- vc[[i]]$g[1]
        P.cov[i, i] <- vc[[i]]$g[1] + attr(vc[[i]], "sc")^2 / nr
      }
      
    }
    
  }
  
  #----------------------------------------------------------------------------
  # Method 1
  # BLUEs and BLUPs
  #----------------------------------------------------------------------------
  
  if (method == 1) {
    # P matrices
    tmp <- dfr.out[, substr(names(dfr.out), 1, 4) == 'blue']
    names(tmp) <- gsub('blue.', '', names(tmp))
    P.cov <- cov(tmp, use = 'pairwise.complete.obs')
    P.cor <- cor(tmp, use = 'pairwise.complete.obs')
    # G matrices
    tmp <- dfr.out[, substr(names(dfr.out), 1, 4) == 'blup']
    names(tmp) <- gsub('blup.', '', names(tmp))
    G.cov <- cov(tmp, use = 'pairwise.complete.obs')
    G.cor <- cor(tmp, use = 'pairwise.complete.obs')
  }

  #----------------------------------------------------------------------------
  # Method 2
  # Based on Z = Y1 + Y2
  #----------------------------------------------------------------------------
  
  if (method == 2) {
    
    if (!is.null(env) & !is.null(rep)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- dfr[, vars[i]] + dfr[, vars[j]]
          model <- lme4::lmer(z ~ (1|g) + (1|g:e) + (1|e/r))
          vcz <- lme4::VarCorr(model) # variance components for z = x + y
          G.cov[i, j] <- G.cov[j, i] <- (vcz$g[1] - G.cov[i, i] - G.cov[j, j]) / 2
          P.cov[i, j] <- P.cov[j, i] <- (vcz$g[1] + vcz$e[1] / ne + attr(vcz, "sc")^2 / ne / nr - P.cov[i, i] - P.cov[j, j]) / 2
        }
      }
    }
    
    if (is.null(env)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- dfr[, vars[i]] + dfr[, vars[j]]
          model <- lme4::lmer(z ~ (1|g) + (1|r))
          vcz <- lme4::VarCorr(model) # variance components for z = x + y
          G.cov[i, j] <- G.cov[j, i] <- (vcz$g[1] - G.cov[i, i] - G.cov[j, j]) / 2
          P.cov[i, j] <- P.cov[j, i] <- (vcz$g[1] + attr(vcz, "sc")^2 / nr - P.cov[i, i] - P.cov[j, j]) / 2
        }
      }
    }
    
    if (is.null(rep)) {
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
          z <- dfr[, vars[i]] + dfr[, vars[j]]
          model <- lme4::lmer(z ~ (1|g) + (1|e))
          vcz <- lme4::VarCorr(model) # variance components for z = x + y
          G.cov[i, j] <- G.cov[j, i] <- (vcz$g[1] - G.cov[i, i] - G.cov[j, j]) / 2
          P.cov[i, j] <- P.cov[j, i] <- (vcz$g[1] + attr(vcz, "sc")^2 / ne - P.cov[i, i] - P.cov[j, j]) / 2
        }
      }
    }
    
    d1 <- diag(diag(G.cov)^{-0.5}, nt, nt)
    d2 <- diag(diag(P.cov)^{-0.5}, nt, nt)
  
    G.cor <- d1 %*% G.cov %*% d1
    P.cor <- d2 %*% P.cov %*% d2
    
  }
  
  #----------------------------------------------------------------------------
  # Method 3
  # Based on average correlation matrix
  #----------------------------------------------------------------------------
  
  if (method == 3) {
    
    if (!is.null(env) & !is.null(rep)) {
      # split by env and rep
      dfr.split <- dfr[, c(vars, env, rep)]
      dfr.split <- split(dfr.split, factor(paste(dfr.split[, env], dfr.split[, rep])))
    }
    
    if (is.null(env)) {
      # split by rep
      dfr.split <- dfr[, c(vars, rep)]
      dfr.split <- split(dfr.split, dfr.split[, rep])
    }
    
    if (is.null(rep)) {
      # split by env
      dfr.split <- dfr[, c(vars, env)]
      dfr.split <- split(dfr.split, dfr.split[, env])
    }
    
    ner <- length(dfr.split) # Number of data frames
    
    cl <- list() # correlation list
    for (i in 1:ner)
      cl[[i]] <- cor(dfr.split[[i]][, 1:nt], use = "pairwise.complete.obs")
    G.cor <- apply(simplify2array(cl), 1:2, mean, na.rm = TRUE)
    P.cor <- cor(dfr[, vars], use = "pairwise.complete.obs")
    
    d1 <- diag(diag(G.cov)^0.5, nt, nt)
    d2 <- diag(diag(P.cov)^0.5, nt, nt)
    
    G.cov <- d1 %*% G.cor %*% d1 # Genotypic covariance matrix
    P.cov <- d2 %*% P.cor %*% d2 # Phenotypic covariance matrix
    
  }
  
  dimnames(G.cov) <- dimnames(P.cov) <- dimnames(G.cor) <- dimnames(P.cor) <- list(vars, vars)

  # results
  
  list(G.cov = G.cov, P.cov = P.cov, G.cor = G.cor, P.cor = P.cor,
       blues = dfr.out[, substr(names(dfr.out), 1, 4) %in% c(geno, 'blue')],
       blups = dfr.out[, substr(names(dfr.out), 1, 4) %in% c(geno, 'blup')])
  
}
