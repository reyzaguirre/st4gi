#' Pesek-Baker Index
#'
#' Function to compute the Pesek-Baker index (Pesek, J. and R.J. Baker., 1969).
#' @param geno The name of the column that identifies the genotypes.
#' @param traits The names of the columns for the traits.
#' @param dgg Desired genetic gains. The default is one standard deviation for each trait.
#' @param units Units for dgg, \code{"sdu"} or \code{"actual"}. See details for more information.
#' @param sf Selected fraction. The default is 0.1.
#' @param dfr The name of the data frame. See details for more information.
#' @param G.cov Genotypic covariance matrix.
#' @param P.cov Phenotypic covariance matrix.
#' @details The Pesek-Baker is an index where relative economic weights have been replaced
#' by desired gains.
#' If \code{dgg} is not specified, the desired genetic gains are set to one
#' standard deviation for each trait. \code{dgg} can be specified in actual units
#' (\code{units = "actual"}) or in standard deviations (\code{units = "sdu"}),
#' defaults to \code{"sdu"}. For example, if you have a trait which is expressed
#' in kilograms and with a standard deviation of 5 kilograms, typing \code{dgg = 2}
#' means a desired genetic gain of 2 standard deviations that corresponds to 10 kilograms.
#' If you type \code{dgg = 2} and \code{units = "actual"} then this means a desired
#' genetic gain of 2 kilograms. If \code{dgg = NULL} then the desired genetic gains
#' will be one standard deviation, no matter if \code{units} is set as \code{"actual"}
#' or \code{"sdu"}.
#' The \code{dfr} must contain the estimations for each genotype with each trait.
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
#' output <- ecm(traits, "geno", "loc", "rep", spg)
#' colnames(output$blups) <- gsub('blup.', '', colnames(output$blups))
#' pesekbaker('geno', traits, dfr = output$blups, G.cov = output$G.cov, P.cov = output$P.cov)
#' # More weight on bc and dm, less on star and nocr
#' dgg <- c(1, 1.5, 1.5, 0.8, 0.8)
#' pesekbaker('geno', traits, dgg, dfr = output$blups, G.cov = output$G.cov, P.cov = output$P.cov)
#' @importFrom stats dnorm qnorm
#' @export

pesekbaker <- function(geno, traits, dgg = NULL, units = c("sdu", "actual"),
                       sf = 0.1, dfr, G.cov, P.cov) {
  
  # Match arguments
  
  units <- match.arg(units)
  
  # Check traits, G.cov, and P.cov are in the same order
  
  if (!identical(traits, rownames(G.cov)) | !identical(traits, rownames(P.cov)))
    warning("Check that traits are in the same order in traits, G.cov and P.cov.")

  # As character
  
  dfr[, geno] <- as.character(dfr[, geno])
  
  # Standard deviations for traits
  
  sdt <- diag(G.cov)^0.5
  
  # Compute index coefficients
  
  if (is.null(dgg)) {
    dgg <- sdt
  } else {
    if (units == "sdu")
      dgg <- dgg * sdt
  }
  
  b <- solve(G.cov) %*% dgg
  dimnames(b) <- list(traits, "coef")
  
  # Response to selection
  
  rs <- NULL
  
  # Compute selection intensity and response to selection
  
  si <- dnorm(qnorm(1 - sf)) / sf
  bPb <- t(b) %*% P.cov %*% b
  for (i in 1:length(traits))
    rs[i] <- si * t(b) %*% G.cov[, i] / sqrt(bPb * G.cov[i, i])
  rsa <- rs * sdt # response to selection in actual units
  names(rs) <- names(rsa)
  
  # Compute index and ranking
  
  m <- as.matrix(dfr[, traits])
  indices <- m %*% b
  dfr <- cbind(dfr, indices)
  colnames(dfr)[colnames(dfr) == 'coef'] <- 'PB.Index'
  dfr$PB.Rank <- rank(-dfr$PB.Index, na.last = "keep")
  dfr <- dfr[order(dfr[, 'PB.Rank']), ]
  
  # Results
  
  list(Desired.Genetic.Gains = dgg,
       Standard.Deviations = sdt,
       Index.Coefficients = b,
       Response.to.Selection = rsa,
       Std.Response.to.Selection = rs,
       Pesek.Baker.Index = dfr)
  
}
