#' Friedman test
#'
#' Performs a Friedman rank sum test with multiple comparisons.
#' @param dfr The name of the data frame.
#' @param y The name of the column for the variable to analyze.
#' @param treat The name of the column that identifies the treatments.
#' @param block The name of the column that identifies the blocks.
#' @param alpha Significant level for comparisons.
#' @param print.mc Logical, if \code{TRUE}, it prints multiple comparisons groups.
#' @param print.pc Logical, if \code{TRUE}, it prints pairwise comparisons tests.
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
#' friedman.t(dfr, "y", "grass", "homeowner")
#' @importFrom stats pchisq
#' @export

friedman.t <- function(dfr, y, treat, block, alpha = 0.05,
                       print.mc = TRUE, print.pc = TRUE) {
  
  # Create a data matrix
  # compute means if more of one evaluation
  
  dmx <- by(dfr[, y], dfr[, c(block, treat)], mean, na.rm = TRUE)
  dmx <- as.data.frame(dmx[, ])
  
  # Check NA's
  
  if (sum(is.na(dmx)) > 0)
    stop("There are NA values")
  
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
 
  # Create groups of non-significant differences
  
  groups <- data.frame(sort(R, decreasing = TRUE), '')
  groups[, 2] <- as.character(groups[, 2])
  
  i <- 1
  l <- 1
  pos <- 1
  
  for (i in 1:(k - 1)) {
    if (pos < k) {
      if (abs(groups[i, 1] - groups[pos + 1, 1]) < lsd) {
        if (i > 1)
          for (m in 1:(i-1))
            groups[m, 2] <- paste0(groups[m, 2], ' ')
        groups[i, 2] <- paste0(groups[i, 2], letters[l])
        for (j in (i + 1):k)
          if (abs(groups[i, 1] - groups[j, 1]) < lsd) {
            groups[j, 2] <- paste0(groups[j, 2], letters[l])
            pos <- j
          } 
        l <- l + 1
      } else {
        if (pos == k - 1 & i == k - 1) {
          for (n in 1:pos)
            groups[n, 2] <- paste0(groups[n, 2], ' ')
          groups[pos + 1, 2] <- paste0(groups[pos + 1, 2], letters[l])
          pos <- pos + 1
        }
      }
    }
  }
  
  colnames(groups) <- c('Sum of ranks', 'Groups')

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
  
  cat("\n     Friedman rank sum test\n\n")
  cat('Chi-square statistic = ', T1, ', df = ', dfn,
      ', p-value = ', p.T1, '\n', sep = '')
  cat('F statistic = ', T2, ', num.df = ', dfn, ', den.df = ', dfd,
      ', p-value = ', p.T2, '\n', sep = '')
  if (print.mc == TRUE & print.pc == TRUE) {
    cat("\n     Multiple comparisons\n\n")
    cat('alpha = ', alpha, ' LSD = ', lsd, '\n\n')
    print(groups)
    cat('\n')
    print(compar)
  }
  if (print.mc == TRUE & print.pc == FALSE) {
    cat("\n     Multiple comparisons\n\n")
    cat('alpha = ', alpha, ' LSD = ', lsd, '\n\n')
    print(groups)
  }
  if (print.mc == FALSE & print.pc == TRUE) {
    cat("\n     Multiple comparisons\n\n")
    cat('alpha = ', alpha, ' LSD = ', lsd, '\n\n')
    print(compar)
  }
  
  output <- list(T1 = T1, p.T1 = p.T1, T2 = T2, p.T2 = p.T2, dfn = dfn, dfd = dfd,
                 lsd = lsd, groups = groups, compar = compar)
  
  invisible(output)
  
}







  
