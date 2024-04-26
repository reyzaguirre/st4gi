#' Pesek-Baker Index
#'
#' Function to compute the Pesek-Baker index (Pesek, J. and R.J. Baker., 1969).
#' @param traits The names of the columns for the traits.
#' @param geno The name of the column that identifies the genotypes.
#' @param env The name of the column that identifies the environments.
#' @param rep The name of the column that identifies the replications.
#' @param dfr The name of the data frame.
#' @param gvals The genotypic values to compute the index, \code{"blups"},
#' \code{"blues"} or \code{"means"}. The default is \code{"blups"}.
#' @param dgg Desired genetic gains. The default is one standard deviation for each trait.
#' @param units Units for dgg, \code{"sdu"} or \code{"actual"}.
#' See details for more information.
#' @param sf Selected fraction. The default is 0.1.
#' @param method The method to compute genotypic covariances. See \code{?ecm} for details.
#' @details The Pesek-Baker is an index where relative economic weights have been replaced
#' by desired gains.
#' 
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
#' @importFrom stats dnorm qnorm
#' @export

pesekbaker <- function(traits, geno, env = NULL, rep, dfr,
                       gvals = c("blups", "blues", "means"),
                       dgg = NULL, units = c("sdu", "actual"),
                       sf = 0.1, method = 1) {

  # Match arguments
  
  gvals <- match.arg(gvals)
  units <- match.arg(units)
  
  # As character
  
  dfr[, geno] <- as.character(dfr[, geno])

  # Number of traits

  nt <- length(traits)
 
  # Run ecm to get covariance and correlation matrices
  
  ccm <- ecm(traits, geno, env, rep, dfr, method)
  
  # Compute gvals
  
  if (gvals == 'blups')
    if (method == 1) {
      dfr.out <- ccm$blups
      } else {
        tmp <- ecm(traits, geno, env, rep, dfr, 1)
        dfr.out <- tmp$blups
      }

  if (gvals == 'blues')
    if (method == 1) {
      dfr.out <- ccm$blues
    } else {
      tmp <- ecm(traits, geno, env, rep, dfr, 1)
      dfr.out <- tmp$blues
    }
  
  if (gvals == "means") {
    dfr.out <- data.frame(geno = unique(dfr[, geno]))
    colnames(dfr.out) <- geno
    tmp <- docomp("mean", traits, c(geno, env), dfr = dfr)
    tmp <- docomp("mean", traits, geno, dfr = tmp)
    dfr.out <- merge(dfr.out, tmp, all = TRUE)
    colnames(dfr.out) <- c("geno", paste("means", traits, sep = "."))
  }
  
  # Compute index
  
  index <- pb_index(traits, dgg, units, sf, nt, dfr.out, ccm$G.cov, ccm$P.cov)
 
  # Results
  
  index
  
}

#' Pesek-Baker Index
#'
#' Function to compute the Pesek-Baker index (Pesek, J. and R.J. Baker., 1969)
#' from BLUEs and BLUPs.
#' @param blues The name of the data frame with BLUEs. See details for more information.
#' @param blups The name of the data frame with BLUPs. See details for more information.
#' @param dgg Desired genetic gains. The default is one standard deviation for each trait.
#' @param units Units for dgg, \code{"sdu"} or \code{"actual"}.
#' See details for more information.
#' @param sf Selected fraction. The default is 0.1.
#' @details The Pesek-Baker is an index where relative economic weights have been replaced
#' by desired gains.
#' The \code{blues} and \code{blups} data frames must have the same structure: 
#' First column for genotypes (with name \code{geno}) and the rest of the columns
#' for the BLUEs and BLUPs for each trait, both data frames with the same order
#' of columns.
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
#' pesekbaker2(output$blues, output$blups)
#' # More weight on bc and dm, less on star and nocr
#' pesekbaker2(output$blues, output$blups, dgg = c(1, 1.5, 1.5, 0.8, 0.8))
#' @export

pesekbaker2 <- function(blues, blups, dgg = NULL,
                        units = c("sdu", "actual"), sf = 0.1) {
  
  # Match arguments
  
  units <- match.arg(units)
  
  # Traits and number of traits
  
  traits <- names(blues)[-1]
  
  nt <- length(traits) 
  
  # Get covariance matrices
  
  P.cov <- cov(blues[, -1], use = 'pairwise.complete.obs')
  G.cov <- cov(blups[, -1], use = 'pairwise.complete.obs')
  
  # Compute index
  
  index <- pb_index(traits, dgg, units, sf, nt, blups, G.cov, P.cov)
  
  # Results
  
  index
  
}

#------------------------------------------------------------------------------
# pb_index computes pesek baker index
#------------------------------------------------------------------------------

pb_index <- function(traits, dgg, units, sf, nt, dfr.out, G.cov, P.cov) {
  
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
  
  si <- dnorm(qnorm(1 - sf)) / sf # selection intensity
  bPb <- t(b) %*% P.cov %*% b
  for (i in 1:nt)
    rs[i] <- si * t(b) %*% G.cov[, i] / sqrt(bPb * G.cov[i, i])
  rsa <- rs * sdt # response to selection in actual units
  names(rs) <- names(rsa)
  
  # Compute index and ranking
  
  m <- as.matrix(dfr.out[, 2:(1 + nt)])
  indices <- m %*% b
  dfr.out <- cbind(dfr.out, indices)
  colnames(dfr.out)[2 + nt] <- "PB.Index"
  dfr.out$PB.Rank <- rank(-dfr.out$PB.Index, na.last = "keep")
  dfr.out <- dfr.out[order(dfr.out[, 'PB.Rank']), ]
  
  # Results
  
  list(Desired.Genetic.Gains = dgg,
       Standard.Deviations = sdt,
       Index.Coefficients = b,
       Response.to.Selection = rsa,
       Std.Response.to.Selection = rs,
       Pesek.Baker.Index = dfr.out)
  
}

