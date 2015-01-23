#'  Pesek-Baker Index
#'
#' Function to compute the Pesek-Baker index (Pesek, J. and R.J. Baker., 1969).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks.
#' @param data The name of the data frame containing the data.
#' @param dgg Desired genetic gains, defaults to one standard deviation.
#' @param units Units for dgg, \code{"actual"} or \code{"sdu"}. See details for more information.
#' @param sf Selected fraction, defaults to 0.1.
#' @author Raul Eyzaguirre
#' @details The Pesek-Baker is an index where relative economic weights have been replaced
#' by desired gains. If \code{dgg} is not specified, the standard deviations of the traits
#' are used. It means that the desired genetic gains are equal to one standard deviation for
#' each trait. \code{dgg} can be specified in actual units (\code{units = "actual"}) or in
#' standard deviations (\code{units = "sdu"}), defaults to \code{"sdu"}. For example,
#' if you have a trait which is expressed in kilograms and with a standard deviation of
#' 5 kilograms, typing \code{dgg = 2} means a desired genetic gain of 2 standard deviations
#' that corresponds to 10 kilograms. If you type \code{dgg = 2} and \code{units = "actual"}
#' then this means a desired genetic gain of 2 kilograms. If \code{dgg = NULL} then the
#' desired genetic gain will be one standard deviation, no matter if \code{units} is set
#' as \code{"actual"} or \code{"sdu"}.
#' To compute the index the package \code{lme4} is needed.
#' @return It returns the desired genetic gains in actual units,
#' the estimated standard deviations,
#' the estimated genetic variances,
#' the estimated correlation matrix,
#' the index coefficients,
#' the response to selection,
#' the standardized response to selection,
#' the Pesek-Baker index value,
#' and the Pesek-Baker index value sorted in descending order.

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
#' # Use 2 standard deviations as desired genetic gains for each trait
#' pesekbaker(c("rytha", "bc", "dm", "star", "nocr"), "geno", "loc", "rep", dgg = 2, spg)
#' @export

pesekbaker <- function(traits, geno, env, rep, data, dgg = NULL, units = "sdu", sf = 0.1) {

  # inits

  gv <- NULL # genetic variance
  pv <- NULL # phenotipic variance
  nt <- length(traits) # number of traits
  ng <- nlevels(factor(data[,geno])) # number of genotypes
  ne <- nlevels(factor(data[,env])) # number of environments
  nr <- nlevels(factor(data[,rep])) # number of replications in each environment
  rs <- NULL # response to selection

  # fitted models by REML

  for (i in 1:nt){
    abc <- data.frame(c1 = data[,traits[i]], c2 = data[,geno], c3 = data[,env], c4 = data[,rep])
    model <- lme4::lmer(c1 ~ (1|c2) + (1|c2:c3) + (1|c3/c4), data = abc)
    gv[i] <- lme4::VarCorr(model)$c2[1]
    pv[i] <- lme4::VarCorr(model)$c2[1] + lme4::VarCorr(model)$'c2:c3'[1]/ne +
      attr(lme4::VarCorr(model), "sc")^2/ne/nr
  }

  # compute correlation and covariance matrices

  df <- data[,c(sapply(traits, c), env, rep)]
  df <- split(df, factor(paste(data[,env], data[,rep]))) # split by env and rep
  corr <- cor(df[[1]][,1:nt], use = "pairwise.complete.obs")
  for (i in 2:length(df))
    corr <- corr + cor(df[[i]][,1:nt], use = "pairwise.complete.obs")
  corr <- corr/length(df)
  S <- diag(gv^.5, nt, nt)
  G <- S%*%corr%*%S
  dimnames(G) <- dimnames(corr)
  P <- G
  diag(P) <- pv

  # compute index coefficients

  if (is.null(dgg) == TRUE) dgg <- gv^.5 else
    if (units == "sdu") dgg <- dgg*gv^.5
  b <- solve(G)%*%dgg
  dimnames(b) <- list(dimnames(corr)[[1]], "coef")

  # response to selection

  si <- dnorm(qnorm(1-sf))/sf # selection intensity
  bPb <- t(b)%*%P%*%b
  for (i in 1:nt)
    rs[i] <- si * t(b)%*%G[,i]/sqrt(bPb*G[i,i])
  rsa <- rs * gv^.5 # response to selection in actual units

  # index calculation

  m <- matrix(NA, ng, nt)
  for (i in 1:nt)
    m[,i] <- tapply(data[,traits[i]], data[,geno], mean, na.rm = T)
  indices <- m %*% b
  rownames(indices) <- levels(data[,geno])
  colnames(indices) <- "PesekBakerIndex"
  orden <- order(indices, decreasing = T)
  sort.ind <- indices
  sort.ind[,1] <- indices[orden]
  rownames(sort.ind) <- levels(data[,geno])[orden]
  colnames(sort.ind) <- "PesekBakerIndex"

  # results

  list(Desired.Genetic.Gains = dgg, Standard.Deviations = gv^.5, Genetic.Variances = gv,
       Correlation.Matrix = corr, Index.Coefficients = b,
       Response.to.Selection = rsa, Std.Response.to.Selection = rs,
       Pesek.Baker.Index = indices, Sorted.Pesek.Baker.Index = sort.ind)
}
