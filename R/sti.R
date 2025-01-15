#' Compute stress tolerance indices
#'
#' Compute several stress tolerance indices.
#' @param dfr The name of the data frame.
#' @param vars The list of variables.
#' @param geno The name of the column that identifies the genotypes.
#' @param normal The identification for the normal group.
#' @param stress The identification for the stress group.
#' @details The names of the columns for the variables must follow the convention
#' \code{y.normal} and \code{y.stress}, so the data frame must be in
#' wide format for the factor with levels \code{normal} and \code{stress}.
#' 
#' For a variable \eqn{y} with values \eqn{y_n} and \eqn{y_s} under normal and
#' stress conditions, the stress intensity (\eqn{si}) is computed:
#' \deqn{
#' si = 1 - \frac{\bar y_s}{\bar y_n}
#' }
#' It ranges between 0 and 1 and the larger the value, the more severe is the stress intensity.
#' Then, the following indices are computed:
#' \itemize{
#'  \item \code{tol} : Tolerance.
#'  \deqn{
#'  tol = y_n - y_s
#'  }  
#'  \item \code{yrr} : Yield reduction ratio.
#'  \deqn{
#'  yrr = 1 - \frac{y_s}{y_n}
#'  }  
#'  \item \code{ssi} : Stress susceptibility index.
#'  \deqn{
#'  ssi = \frac{yrr}{si}
#'  }  
#'  \item \code{sti} : Stress tolerance index.
#'  \deqn{
#'  sti = \frac{y_n \times y_s}{\bar y_n^2}
#'  }
#'  \item \code{mp}  : Mean productivity.
#'  \deqn{
#'  mp = \frac{y_n + y_s}{2}
#'  }  
#'  \item \code{gmp} : Geometric mean productivity.
#'  \deqn{
#'  gmp = \sqrt{y_n \times y_s}
#'  }  
#'  }
#' @return It returns a data frame with the indices and a data frame with
#' the stress intensity values for each variable. The names for the indices
#' follow the convention \code{y.index}.
#' @author Raul Eyzaguirre.
#' @references
#' Fernandez, G.C.J. (1992). Effective Selection Criteria for Assessing Stress Tolerance.
#' In: Kuo, C.G., Ed., Proceedings of the International Symposium on Adaptation of
#' Vegetables and Other Food Crops in Temperature and Water Stress, AVRDC Publication,
#' Tainan, 257-270.
#' @examples
#' vars <- c("nmtp", "mtwp", "nnomtp")
#' sti(potatostress, vars, 'genotype', 'DTWW', 'DTWS')
#' @export

sti <- function(dfr, vars, geno, normal, stress) {
  
  # Variable names
  
  vars.normal <- paste(vars, normal, sep = '.')
  vars.stress <- paste(vars, stress, sep = '.')
  
  # Stress intensity
  
  si.out <- NULL
  si.name <- NULL
  
  # Indices 
  
  for (i in 1:length(vars)) {
    
    # Normal value of the variable

    y_n <- dfr[, vars.normal[i]]
    
    # Stress value of the variable
    
    y_s <- dfr[, vars.stress[i]]
    
    # Tolerance
    
    index.name <- paste(vars[i], 'tol', sep = '.')
    dfr[, index.name] <- y_n - y_s
    
    # Yield reduction ratio
    
    index.name <- paste(vars[i], 'yrr', sep = '.')
    dfr[, index.name] <- 1 - y_s / y_n
    tmp <- index.name # Keep name for next index
    
    # Stress intensity
    
    si <- 1 - mean(y_s) / mean(y_n)
    si.out <- c(si.out, si)
    si.name <- c(si.name, vars[i])
    
    # Stress susceptibility index
    
    index.name <- paste(vars[i], 'ssi', sep = '.')
    dfr[, index.name] <- dfr[, tmp] / si
    
    # Stress tolerance index
    
    index.name <- paste(vars[i], 'sti', sep = '.')
    dfr[, index.name] <- y_n * y_s / mean(y_n)^2

    # Mean productivity
    
    index.name <- paste(vars[i], 'mp', sep = '.')
    dfr[, index.name] <- (y_n + y_s) / 2

    # Geometric mean productivity
    
    index.name <- paste(vars[i], 'gmp', sep = '.')
    dfr[, index.name] <- sqrt(y_n * y_s)
   
  }
  
  # Return
  
  si.values <- data.frame(variable = si.name,
                          si = si.out)
    
  list(index.dfr = dfr,
       si.values = si.values)
  
}
