#'  Pesek-Baker Index
#'
#' Function to compute the Pesek-Baker index (Pesek, J. and R.J. Baker., 1969).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @param means The genotypic means to compute the index, \code{"single"}
#' or \code{"fitted"}. The default is \code{"single"}. See details for more information.
#' @param dgg Desired genetic gains. The default is one standard deviation for each trait.
#' @param units Units for dgg, \code{"actual"} or \code{"sdu"}. See details for more information.
#' @param sf Selected fraction. The default is 0.1.
#' @param method The method to compute genotypic covariances. See ?ecm for details.
#' @details The Pesek-Baker is an index where relative economic weights have been replaced
#' by desired gains.
#' 
#' If \code{means = "fitted"} then fitted means under the model are used for each genotype.
#' Otherwise single arithmetic means are computed over the replications for each genotype
#' at each environment and then for each genotype over environments.
#' 
#' If \code{dgg} is not specified, the genotypic standard deviations of the traits are used.
#' It means that the desired genetic gains are equal to one standard deviation for
#' each trait. \code{dgg} can be specified in actual units (\code{units = "actual"}) or in
#' standard deviations (\code{units = "sdu"}), defaults to \code{"sdu"}. For example,
#' if you have a trait which is expressed in kilograms and with a standard deviation of
#' 5 kilograms, typing \code{dgg = 2} means a desired genetic gain of 2 standard deviations
#' that corresponds to 10 kilograms. If you type \code{dgg = 2} and \code{units = "actual"}
#' then this means a desired genetic gain of 2 kilograms. If \code{dgg = NULL} then the
#' desired genetic gains will be one standard deviation, no matter if \code{units} is set
#' as \code{"actual"} or \code{"sdu"}.
#' @return It returns:
#' \itemize{
#' \item \code{$Desired.Genetic.Gains}, the desired genetic gains in actual units,
#' \item \code{$Standard.Deviations}, the estimated genotypic standard deviations,
#' \item \code{$Index.Coefficients}, the index coefficients,
#' \item \code{$Response.to.Selection}, the response to selection,
#' \item \code{$Std.Response.to.Selection}, the standardized response to selection, and
#' \item \code{$Pesek.Baker.Index}, a data frame with the genotypic means for each trait,
#' the Pesek-Baker index, and the rank for each genotype according to the index.
#' }
#' @author Raul Eyzaguirre.
#' @references
#' Pesek, J. and R.J. Baker.(1969). Desired improvement in relation to selection indices.
#' Can. J. Plant. Sci. 9:803-804.
#' @examples
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' pesekbaker(traits, "geno", "loc", "rep", spg)
#' # More weight on bc and dm, less on star and nocr
#' pesekbaker(traits, "geno", "loc", "rep", spg, dgg = c(1, 1.5, 1.5, 0.8, 0.8))
#' @importFrom stats cor dnorm qnorm
#' @export

pesekbaker <- function(traits, geno, env = NULL, rep, dfr, means = c("single", "fitted"),
                       dgg = NULL, units = c("sdu", "actual"), sf = 0.1, method = 1) {

  # match arguments
  
  means <- match.arg(means)
  units <- match.arg(units)
  
  # As character
  
  dfr[, geno] <- as.character(dfr[, geno])

  # inits

  nt <- length(traits) # number of traits
  rs <- NULL # response to selection
  
  # run ecm to get covariance and correlationi matrices
  
  cm <- ecm(traits, geno, env, rep, dfr, method)
  
  # standard deviations for traits
  
  sdt <- diag(cm$G.Cov)^0.5

  # compute index coefficients

  if (is.null(dgg)) {
    dgg <- sdt
  } else {
    if (units == "sdu")
      dgg <- dgg * sdt
  }
  
  b <- solve(cm$G.Cov) %*% dgg
  dimnames(b) <- list(traits, "coef")

  # response to selection

  si <- dnorm(qnorm(1 - sf)) / sf # selection intensity
  bPb <- t(b) %*% cm$P.Cov %*% b
  for (i in 1:nt)
    rs[i] <- si * t(b) %*% cm$G.Cov[, i] / sqrt(bPb * cm$G.Cov[i, i])
  rsa <- rs * sdt # response to selection in actual units
  names(rs) <- names(rsa)

  # index calculation
  
  dfr.out <- data.frame(geno = unique(dfr[, geno]))
  colnames(dfr.out) <- geno
  
  if (means == "single") {
    temp <- docomp("mean", traits, c(geno, env), dfr = dfr)
    temp <- docomp("mean", traits, geno, dfr = temp)
    dfr.out <- merge(dfr.out, temp, all = TRUE)
    colnames(dfr.out) <- c("geno", paste("m", traits, sep = "."))
  }
  
  if (means == "fitted") {
    for (i in 1:nt) {
      if (!is.null(env)) {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", geno, ":", env,
                               ") + (1|", env, "/", rep, ")"))
        fm <- lme4::lmer(ff, dfr)
        }
      if (is.null(env)) {
        ff <- as.formula(paste(traits[i], "~", geno, "- 1 + (1|", rep, ")"))
        fm <- lme4::lmer(ff, dfr)
      }
      temp <- as.data.frame(lme4::fixef(fm))
      colnames(temp) <- paste("f", traits[i], sep = ".")
      temp[, geno] <- substring(rownames(temp), nchar(geno) + 1)
      dfr.out <- merge(dfr.out, temp, all = TRUE)
    }
  }

  m <- as.matrix(dfr.out[, 2:(1 + nt)])
  indices <- m %*% b
  dfr.out <- cbind(dfr.out, indices)
  colnames(dfr.out)[2 + nt] <- "PB.Index"
  dfr.out$PB.Rank <- rank(-dfr.out$PB.Index, na.last = "keep")
  
  # results

  list(Desired.Genetic.Gains = dgg,
       Standard.Deviations = sdt,
       Index.Coefficients = b,
       Response.to.Selection = rsa,
       Std.Response.to.Selection = rs,
       Pesek.Baker.Index = dfr.out)
}
