# #' Adjust values following the method of Westcott
# #'
# #' This function adjust the observed values of an experiment planted following
# #' the method described by Westcott (1981) with a grid of checks.
# #' @param y The name of the column for the variable to adjust.
# #' @param geno The name of the column that identifies the genotypes.
# #' @param ck1 Name of check 1.
# #' @param ck2 Name of check 2.
# #' @param row The name of the column that identifies the rows.
# #' @param col The name of the column that identifies the columns.
# #' @param ncb Number of columns between two check columns.
# #' @param nrs Number of rows to span the row of the plot.
# #' @param method The method to fit the values. See details.
# #' @param ind Logical. See details.
# #' @param p The proportion of the check values differences used for the adjustment.
# #' See details.
# #' @param dfr The name of the data frame.
# #' @details The values of the selected variable \code{y} are adjusted using some mean
# #' of the values of all the checks located on the row of the plot plus the \code{nrs}
# #' rows at each side of the row of the plot. If \code{method = "flat"} the simple mean
# #' of the checks is used for the adjustment, if \code{method = "weighted"} a weighted mean
# #' is used for the adjustmen, where the checks closer to the plot get more weight.
# #' If \code{p = 1} then the values are adjusted in the same proportion that the checks
# #' vary around the field, for values lower than 1 the values are adjusted based on that
# #' proportion over the checks variation, if \code{p = 0} then there is no adjustment.
# #' If \code{ind = TRUE}, each check is centered around its own mean before the adjustment,
# #' if \code{ind = FALSE}, then both checks are centered around the overall mean.
# #' If \code{method = "flat"}, \code{ncb = 5}, \code{nrs = 1}, and \code{ind = FALSE}, then
# #' it corresponds to the method proposed by Westcott.
# #' If not specified, \code{nrs = floor(ncb / 2)}.
# #' If the layout does not correspond with the Westcott method, then the observed values
# #' are adjusted with the values of the checks planted nearby (this is in the rectangular
# #' region definied by \code{floor(ncb / 2)} columns and rows at each side of
# #' the plot) and a warning is issued.
# #' @return It returns the adjusted values.
# #' @author Raul Eyzaguirre.
# #' @references
# #' Westcott, B. (1981). Two methods for early generation yield assessment in winter wheat.
# #' In: Proc. of the 4th meeting of the Biometrics in Plant Breeding Section of Eucarpia.
# #' INRA Poitier, France, pp 91-95.
# #' @examples 
# #' # Create design
# #' dfr <- cr.w(1:1000, "Dag", "Cem", 40, 10)
# #' dfr <- dfr$book
# #' # Create some random data
# #' dfr$y <- rnorm(dim(dfr)[1])
# #' # Run adjustment
# #' aj.w("y", "geno", "Dag", "Cem", "row", "col", dfr = dfr)
# #' @export
# aj.w <- function(y, geno, ck1, ck2, row, col, ncb = 10, nrs = NULL,
#                  method = c("weighted", "flat"), ind = TRUE, p = 0.5, dfr) {
#   
#   # Match arguments
#   
#   method <- match.arg(method)
#   
#   # Define nrs
#   
#   if (is.null(nrs))
#     nrs <- floor(ncb / 2)
#   
#   # Error and warning messages
#   
#   out <- ck.pos(row, col, NULL, dfr)
#   
#   if (out$nplot > 0)
#     stop("More than one genotype in the same position. Run check.pos to look over.")
#   
#   out <- ck.w(y, geno, ck1, ck2, row, col, ncb, dfr)
#   eval.cond <- sum(out$c1, out$c2, out$c3, out$c4, out$c5)
#   
#   if (out$c1 == 1)
#     warning("There are plots in the columns of checks with other genotypes planted.")
#   
#   if (out$c2 == 1)
#     warning("There are plots in the columns of genotypes with checks planted.")
#   
#   if (out$c3 == 1)
#     warning("There are columns of checks with missing plots.")
# 
#   if (out$c4 == 1)
#     warning("There are columns of checks without alternating checks.")
#   
#   if (out$c5 == 1)
#     warning("There are plots with genotypes without a check plot to the left or to the right.")
#   
#   if (eval.cond > 0)
#     warning("Adjusted values are obtained with the values of the checks nearby.")
# 
#   # Get a copy of y for the adjusted values
#   
#   y.w <- paste(y, 'w', sep = '.') 
#   dfr[, y.w] <- dfr[, y]
# 
#   # Compute means for checks
#   
#   if (ind) {
#     ck1.mean <- mean(dfr[dfr[, geno] == ck1, y], na.rm = TRUE)
#     ck2.mean <- mean(dfr[dfr[, geno] == ck2, y], na.rm = TRUE)
#   } else {
#     ck1.mean <- mean(dfr[dfr[, geno] %in% c(ck1, ck2), y], na.rm = TRUE)
#     ck2.mean <- mean(dfr[dfr[, geno] %in% c(ck1, ck2), y], na.rm = TRUE)
#   }
#   
#   # Center check values and compute relative variation adjusted by p
#   
#   cond1 <- dfr[, geno] == ck1
#   cond2 <- dfr[, geno] == ck2
# 
#   dfr[cond1, y.w] <- (dfr[cond1, y.w] - ck1.mean) / ck1.mean * p
#   dfr[cond2, y.w] <- (dfr[cond2, y.w] - ck2.mean) / ck2.mean * p
#   
#   # Replace missing values with 0 (this is the centered mean)
#   
#   dfr[dfr[, geno] %in% c(ck1, ck2) & is.na(dfr[, y.w]), y.w] <- 0
#   
#   # Choose function for adjustment
#   
#   if (eval.cond == 0 & method == "weighted")
#     foo <- foo.weig
#   if (eval.cond == 0 & method == "flat")
#     foo <- foo.flat
#   if (eval.cond > 0) {
#     nrs <- floor(ncb / 2)
#     ncb <- floor(ncb / 2)
#     foo <- foo.flat
#   }
#   
#   # Make adjustment
#   
#   for (i in 1:dim(dfr)[1]) {
#     if (dfr[i, geno] %in% c(ck1, ck2)) {
#       af <- 0
#     } else {
#         af <- foo(i, y.w, geno, ck1, ck2, row, col, ncb, nrs, dfr)
#     }
#     dfr[i, y.w] <- dfr[i, y.w] / (1 + af)
#   }
#   
#   # Return
#   
#   dfr
#   
# }
# 
# # A function for Westcott adjustment with method 2
# 
# foo.weig <- function(x, y.w, geno, ck1, ck2, row, col, ncb, nrs, dfr) {
#   
#   # Identify row and column for plot
#   
#   geno.row <- dfr[x, row]
#   geno.col <- dfr[x, col]
#   
#   # Identify columns for checks (left and right)
#   
#   columns <- (geno.col - ncb):(geno.col + ncb)
#   tmp <- dfr[dfr[, row] == geno.row & dfr[, col] %in% columns & dfr[, geno] %in% c(ck1, ck2), ]
#   col.lf <- min(tmp[, col])
#   col.rg <- max(tmp[, col])
#   
#   # Identify rows for checks and define weights
#   
#   row.ck <- (geno.row - nrs):(geno.row + nrs)
#   row.wg <- c(1:nrs, nrs + 1, nrs:1)
#   
#   # Delete nonexistent rows
#   
#   tmp <- dfr[dfr[, col] == geno.col, ]
#   valid.rows <- row.ck %in% tmp[, row]
#   row.ck <- row.ck[valid.rows]
#   row.wg <- row.wg[valid.rows]
#   
#   # Check values on the left and right
#   
#   ck.lf <- dfr[dfr[, row] %in% row.ck & dfr[, col] == col.lf, c(row, y.w)]
#   ck.rg <- dfr[dfr[, row] %in% row.ck & dfr[, col] == col.rg, c(row, y.w)]
#   
#   # Sort by row
#   
#   ck.lf <- ck.lf[sort.int(ck.lf[, row], index.return = TRUE)$ix, ]
#   ck.rg <- ck.rg[sort.int(ck.rg[, row], index.return = TRUE)$ix, ]
#   
#   # Get left and right means for adjustment
#   
#   m.ck.lf <- sum(ck.lf[, y.w] * row.wg) / sum(row.wg)
#   m.ck.rg <- sum(ck.rg[, y.w] * row.wg) / sum(row.wg)
#   
#   # Get adjustment factor
#   
#   af <- m.ck.lf + (geno.col - col.lf) * (m.ck.rg - m.ck.lf) / (col.rg - col.lf)
#   
#   # Return
#   
#   af
#   
# }
# 
# # A function for Westcott adjustment with method 1 or
# # when layout does not follow the Westcott method
# 
# foo.flat <- function(x, y.w, geno, ck1, ck2, row, col, ncb, nrs, dfr) {
#   
#   # Identify row and column for plot
#   
#   geno.row <- dfr[x, row]
#   geno.col <- dfr[x, col]
#   
#   # Identify all check values in the neighbourhood
#   
#   rows <- (geno.row - nrs):(geno.row + nrs)
#   columns <- (geno.col - ncb):(geno.col + ncb)
# 
#   tmp <- dfr[dfr[, row] %in% rows & dfr[, col] %in% columns & dfr[, geno] %in% c(ck1, ck2), ]
# 
#   # Get adjustment factor
#   
#   af <- mean(tmp[, y.w])
#   
#   # Return
#   
#   af
#   
# }

# -----------------------------------------------------------------------------

# #' Optimize Westcott adjustment
# #'
# #' This function optimizes the parameter values of the Westcott adjustment.
# #' @param y The name of the column for the variable to adjust.
# #' @param geno The name of the column that identifies the genotypes.
# #' @param replicated A \code{yes/no} column to identify replicated genotypes.
# #' @param ck1 Name of check 1.
# #' @param ck2 Name of check 2.
# #' @param row The name of the column that identifies the rows.
# #' @param col The name of the column that identifies the columns.
# #' @param ncb Number of columns between two check columns.
# #' @param nrs Set of values for the \code{nrs} param on \code{aj.w}.
# #' @param method Set of values for the \code{method} param on \code{aj.w}.
# #' @param ind Set of values for the \code{ind} param on \code{aj.w}.
# #' @param p Set of values for the \code{p} param on \code{aj.w}.
# #' @param opt Criteria for optimization, \code{F.value}, \code{MSE} or \code{CV}.
# #' See details.
# #' @param dfr The name of the data frame.
# #' @details If \code{optim = F.value} then the F value for ANOVA is maximized.
# #' If \code{optim = MSE} then the MSE is minimized.
# #' If \code{optim = CV} then the CV is minimized.
# #' @return It returns optimal set of values under the selected criteria.
# #' @author Raul Eyzaguirre.
# #' @export
# 
# optim.w <- function(y, geno, replicated, ck1, ck2, row, col, ncb, nrs, method, ind,
#                     p, opt = c("F.value", "MSE", "CV"), dfr) {
#   
#   # Match arguments
#   
#   opt <- match.arg(opt)
# 
#   # All settings to run
#   
#   nrs.l <- length(nrs)
#   method.l <- length(method)
#   ind.l <- length(ind)
#   p.l <- length(p)
#   
#   total <- nrs.l * method.l * ind.l * p.l
#   
#   # Anova without adjustment
#   
#   tmp <- dfr[dfr[, replicated] == 'yes', ]
#   ff <- as.formula(paste(y, '~', geno))
#   at <- anova(aov(ff, tmp))
#   MSE <- at[2, 3]
#   F.value <- at[1, 4]
#   CV <- sqrt(at[2, 3]) / mean(tmp[, y], na.rm = T) * 100
#   
#   # Original values
#   
#   best.method <- ''
#   best.ind <- ''
#   best.nrs <- ''
#   best.p <- ''
#   best.MSE <- MSE
#   best.F <- F.value 
#   best.CV <- CV
#   
#   # Run adjustments
#   
#   y.aj <- paste0(y, '.aj')
#   i <- 0
#   
#   for (method.v in method[1:method.l])
#     for (ind.v in ind[1:ind.l])
#       for (nrs.v in nrs[1:nrs.l])
#         for (p.v in p[1:p.l]) {
# 
#           dfr <- aj.w(y, geno, ck1, ck2, row, col, ncb, nrs.v, method.v, ind.v, p.v, dfr)
#           tmp <- dfr[dfr[, replicated] == 'yes', ]
#           ff <- as.formula(paste(y.aj, '~', geno))
#           at.aj <- anova(aov(ff, tmp))
#           
#           new.F <- at.aj[1, 4]
#           new.MSE <- at.aj[2, 3]
#           new.CV <- sqrt(at.aj[2, 3]) / mean(tmp[, y.aj], na.rm = T) * 100
#           
#           cond1 <- opt == "F.value" & new.F > best.F
#           cond2 <- opt == "MSE" & new.MSE < best.MSE
#           cond3 <- opt == "CV" & new.CV < best.CV
#           
#           if (cond1 | cond2 | cond3) {
#             best.method <- method.v
#             best.ind <- ind.v
#             best.nrs <- nrs.v
#             best.p <- p.v
#             best.MSE <- new.MSE
#             best.F <- new.F
#             best.CV <- new.CV
#             cat('\n-------------------------------\n')
#             cat("Method = ", best.method, '\n')
#             cat("Ind    = ", best.ind, '\n')
#             cat("nrs    = ", best.nrs, '\n')
#             cat("p      = ", best.p, '\n')
#             cat("MSE    = ", MSE, '\n')
#             cat("MSE.aj = ", best.MSE, '\n')
#             cat("F      = ", F.value, '\n')
#             cat("F.aj   = ", best.F, '\n')
#             cat("CV     = ", CV, '\n')
#             cat("CV.aj  = ", best.CV, '\n')
#             cat('-------------------------------\n')
#           }
#           
#           i <- i + 1
#           cat(round(i / total * 100, 0), '% ', sep = "")
#           
#         }
#   
# }




