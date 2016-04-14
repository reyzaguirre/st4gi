#'  Pesek-Baker Index
#'
#' Function to compute the Pesek-Baker index (Pesek, J. and R.J. Baker., 1969).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications. Must be defined if \code{model = "gxe"}.
#' @param data The name of the data frame containing the data.
#' @param means The genotypic means to compute the index, \code{"single"}
#' or \code{"fitted"}. The default is \code{"single"}. See details for more information.
#' @param model Type of model, \code{"gxe"} for a model with gxe interaction or \code{"g+e"}
#' for a model without interaction. The default is \code{"gxe"}. See details for more information.
#' @param dgg Desired genetic gains. The default is one standard deviation for each trait.
#' @param units Units for dgg, \code{"actual"} or \code{"sdu"}. See details for more information.
#' @param sf Selected fraction. The default is 0.1.
#' @author Raul Eyzaguirre
#' @details The Pesek-Baker is an index where relative economic weights have been replaced
#' by desired gains.
#' 
#' By default a model with components for genotypes, environments, genotypes by environments
#' interaction and replications nested into environments is fitted (\code{model = "gxe"}).
#' If \code{model = "g+e"} then a model with components for genotypes
#' and environments is fitted, and in this case the gxe variance includes the gxe plus the
#' error variance. Response to selection is only computed when \code{model = "gxe"}.
#' 
#' If \code{means = "fitted"} then the model specified in \code{model} is used to fit
#' the means of the genotypes. Otherwise single arithmetic means are computed over the
#' replications for each genotype at each environment and then for each genotype over environments.
#' 
#' If \code{dgg} is not specified, the standard deviations of the traits are used.
#' It means that the desired genetic gains are equal to one standard deviation for
#' each trait. \code{dgg} can be specified in actual units (\code{units = "actual"}) or in
#' standard deviations (\code{units = "sdu"}), defaults to \code{"sdu"}. For example,
#' if you have a trait which is expressed in kilograms and with a standard deviation of
#' 5 kilograms, typing \code{dgg = 2} means a desired genetic gain of 2 standard deviations
#' that corresponds to 10 kilograms. If you type \code{dgg = 2} and \code{units = "actual"}
#' then this means a desired genetic gain of 2 kilograms. If \code{dgg = NULL} then the
#' desired genetic gain will be one standard deviation, no matter if \code{units} is set
#' as \code{"actual"} or \code{"sdu"}.
#' 
#' @return It returns:
#' \itemize{
#' \item \code{$Desired.Genetic.Gains}, the desired genetic gains in actual units,
#' \item \code{$Standard.Deviations}, the estimated standard deviations,
#' \item \code{$Genetic.Variances}, the estimated genetic variances,
#' \item \code{$Correlation.Matrix}, the estimated correlation matrix,
#' \item \code{$Index.Coefficients}, the index coefficients,
#' \item \code{$Response.to.Selection}, the response to selection,
#' \item \code{$Std.Response.to.Selection}, the standardized response to selection, and
#' \item \code{$Pesek.Baker.Index}, a data frame with the genotypic means for each trait,
#' the Pesek-Baker index, and the rank for each genotype according to the index.
#' }
#' @references
#' Pesek, J. and R.J. Baker.(1969). Desired improvement in relation to selection indices.
#' Can. J. Plant. Sci. 9:803-804.
#' @examples
#' # The data
#' head(spg)
#' str(spg)
#'
#' # Run Pesek-Baker index with all the traits
#' pesekbaker(c("rytha", "bc", "dm", "star", "nocr"), "geno", "loc", "rep", spg)
#'
#' # Use different desired genetic gains for each trait,
#' # more weight on bc and dm, less on star and nocr.
#' pesekbaker(c("rytha", "bc", "dm", "star", "nocr"), "geno", "loc", "rep", spg,
#'            dgg = c(1, 1.5, 1.5, 0.8, 0.8))
#' @export

pesekbaker <- function(traits, geno, env, rep = NULL, data, means = "single",
                       model = "gxe", dgg = NULL, units = "sdu", sf = 0.1) {

  # inits

  gv <- NULL # genetic variance
  pv <- NULL # phenotipic variance
  nt <- length(traits) # number of traits
  ng <- nlevels(factor(data[, geno])) # number of genotypes
  ne <- nlevels(factor(data[, env])) # number of environments
  if (!is.null(rep))
    nr <- nlevels(factor(data[, rep])) # number of replications in each environment
  rs <- NULL # response to selection

  # fitted models by REML for variance components

  if (model == "gxe") {
    for (i in 1:nt) {
      ff <- as.formula(paste(traits[i], "~ (1|", geno, ") + (1|", geno, ":", env,
                             ") + (1|", env, "/", rep, ")"))
      fm <- lme4::lmer(ff, data = data)
      gv[i] <- lme4::VarCorr(fm)[[2]][1]
      pv[i] <- lme4::VarCorr(fm)[[2]][1] + lme4::VarCorr(fm)[[1]][1] / ne +
        attr(lme4::VarCorr(fm), "sc")^2 / ne / nr
    }
  }
  if (model == "g+e") {
    for (i in 1:nt) {
      ff <- as.formula(paste(traits[i], "~ (1|", geno, ") + (1|", env, ")"))
      fm <- lme4::lmer(ff, data = data)
      gv[i] <- lme4::VarCorr(fm)[[1]][1]
    }
  }

  # compute correlation and covariance matrices

  if (!is.null(rep)) {
    df <- data[, c(sapply(traits, c), env, rep)]
    df <- split(df, factor(paste(data[, env], data[, rep]))) # split by env and rep
  }
  if (is.null(rep)) {
    df <- data[, c(sapply(traits, c), env)]
    df <- split(df, data[, env]) # split by env
  }
  ner <- length(df)
  ll <- paste("cor", 1:ner, sep = "_")
  my.list <- list()
  for (i in 1:ner)
    my.list[[ll[i]]] <- cor(df[[i]][, 1:nt], use = "pairwise.complete.obs")
  corr <- apply(simplify2array(my.list), 1:2, mean, na.rm = TRUE)
  
  S <- diag(gv^0.5, nt, nt)
  G <- S %*% corr %*% S
  dimnames(G) <- dimnames(corr)
  if (model == "gxe") {
    P <- G
    diag(P) <- pv
  }

  # compute index coefficients

  if (is.null(dgg)) {
    dgg <- gv^0.5
  } else {
      if (units == "sdu") dgg <- dgg * gv^0.5
  }
  b <- solve(G) %*% dgg
  dimnames(b) <- list(dimnames(corr)[[1]], "coef")

  # response to selection

  if (model == "gxe") {
    si <- dnorm(qnorm(1 - sf)) / sf # selection intensity
    bPb <- t(b) %*% P %*% b
    for (i in 1:nt)
      rs[i] <- si * t(b) %*% G[, i] / sqrt(bPb * G[i, i])
    rsa <- rs * gv^0.5 # response to selection in actual units
  } else {
    rsa <- "NA"
    rs <- "NA"
  }

  # index calculation
  
  outind <- data.frame(geno = levels(factor(data[, geno])))
  colnames(outind) <- geno
  
  if (means == "single") {
    temp <- docomp("mean", traits, c(geno, env), data = data)
    temp <- docomp("mean", traits, geno, data = temp)
    outind <- merge(outind, temp, all = TRUE)
    colnames(outind) <- c("geno", paste("m", traits, sep = "."))
  }
  
  if (means == "fitted") {
    for (i in 1:nt) {
      if (model == "gxe") {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", geno, ":", env,
                               ") + (1|", env, "/", rep, ")"))
        fm <- lme4::lmer(ff, data = data)
        }
      if (model == "g+e") {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", env, ")"))
        fm <- lme4::lmer(ff, data = data)
      }
      temp <- as.data.frame(lme4::fixef(fm))
      colnames(temp) <- paste("f", traits[i], sep = ".")
      temp$geno <- substring(rownames(temp), 5)
      outind <- merge(outind, temp, all = TRUE)
    }
  }

  m <- as.matrix(outind[, 2:(1 + nt)])
  indices <- m %*% b
  outind <- cbind(outind, indices)
  colnames(outind)[2 + nt] <- "PB.Index"
  outind$PB.Rank <- rank(-outind$PB.Index, na.last = "keep")
  
  # results

  list(Desired.Genetic.Gains = dgg, Standard.Deviations = gv^0.5, Genetic.Variances = gv,
       Correlation.Matrix = corr, Index.Coefficients = b,
       Response.to.Selection = rsa, Std.Response.to.Selection = rs,
       Pesek.Baker.Index = outind)
}
