#' Friedman test
#'
#' Performs a Friedman rank sum test with multiple comparisons.
#' @param trait The trait to analyze.
#' @param treat The treatments.
#' @param block The blocks.
#' @param dfr The name of the data frame.
#' @param alpha Significant level for comparisons.
#' @param print.text Logical, if \code{TRUE}, it prints output text to console.
#' @details It tests the null hypothesis that the treatments have identical effects.
#' @return It returns the Friedman test and a the multiple comparisons.
#' @author Raul Eyzaguirre.
#' @references
#' Conover, W. J.(1999). Practical Nonparametric Statistics..
#' @examples
#' # From Conover (1999)
#' dfr <- data.frame(y = c(4, 4, 3, 3, 4, 2, 1, 2, 3.5, 4, 4, 3.5,
#'                         3, 2, 1.5, 1, 2, 2, 3, 4, 1, 1, 2, 1,
#'                         2, 3.0, 1.5, 2, 1, 2, 2, 1, 2, 3, 3, 2,
#'                         1, 1, 4, 4, 3, 4, 4, 3, 3.5, 2, 1, 3.5),
#'                   grass = gl(4, 12),
#'                   homeowner = rep(1:12, 4))
#' friedman.t("y", "grass", "homeowner", dfr)
#' @importFrom stats pchisq
#' @export

friedman.t <- function(trait, treat, block, dfr, alpha = 0.05, print.text = TRUE) {
  
  # Create a data matrix
  # compute means if more of one evaluation
  
  dmx <- by(dfr[, trait], dfr[, c(block, treat)], mean, na.rm = TRUE)
  dmx <- as.data.frame(dmx[, ])
  
  # Number of treatments and blocks
  
  b <- dim(dmx)[1]
  k <- dim(dmx)[2]
  
  # Degrees of freedom
  
  dfn <- k - 1
  dfd <- (b - 1) * (k - 1)
  
  # Compute ranks matrix
  
  rmx <- t(apply(dmx, 1, rank))
  
  # Sum of ranks for each treatment
  
  R <- apply(rmx, 2, sum)
  
  # Compute T1 statistic (adjusted by ties)
  
  A1 <- sum(rmx^2)
  C1 <- b * k * (k + 1)^2 / 4
  T1 <- dfn * (sum(R^2) - b * C1) / (A1 - C1)
  
  # Compute T2 statistic
  
  T2 <- (b - 1) * T1 / (b * dfn - T1)
  
  # p-values for T1 and T2
  
  p.T1 <- 1 - pchisq(T1, dfn)
  p.T2 <- 1 - pf(T2, dfn, dfd)
  
  # Compute LSD value with alpha value
  
  std.diff <- (2 * (b * A1 - sum(R^2)) / dfd)^0.5
  lsd <- qt(1 - alpha/2, dfd) * std.diff
  
  # Generate all pairwise comparisons
  
  comb <- utils::combn(k, 2)
  diff <- R[comb[1, ]] - R[comb[2, ]]
  t.val <- diff / std.diff
  p.t <- (1 - pt(abs(t.val), dfd)) * 2
  lcl <- diff - lsd
  ucl <- diff + lsd
  
  # A data.frame for comparisons
  
  tnames <- colnames(dmx)
  
  compar <- data.frame(difference = diff,
                       pvalue = round(p.t, 4),
                       signif = NA,
                       LCL = round(lcl, 2),
                       UPL = round(ucl, 2))
  
  # Add labels for significance
  
  compar$signif <- "ns"
  compar$signif[compar$pvalue <= 0.1] <- "."
  compar$signif[compar$pvalue <= 0.05] <- "*"
  compar$signif[compar$pvalue <= 0.01] <- "**"
  compar$signif[compar$pvalue <= 0.001] <- "***"
  
  # Informative row names for compar
  
  rownames(compar) <- paste(tnames[comb[1, ]], '-', tnames[comb[2, ]])
  
  # Output
  
  if (print.text == TRUE) {
    cat("\n     Friedman rank sum test\n\n")
    cat('Chi-square statistic = ', T1, ', df = ', dfn,
        ', p-value = ', p.T1, '\n\n', sep = '')
    cat('F statistic = ', T2, ', num.df = ', dfn, ', den.df = ', dfd,
        ', p-value = ', p.T2, '\n\n', sep = '')
    cat("\n     Multiple comparisons\n\n")
    cat('alpha = ', alpha, '\n')
    print(compar)
  }

  output <- list(T1 = T1, p.T1 = p.T1, T2 = T2, p.T2 = p.T2, dfn = dfn, dfd = dfd,
                 lsd = lsd, compar = compar)
  
  invisible(output)
  
}







  
