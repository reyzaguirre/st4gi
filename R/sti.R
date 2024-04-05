#' Compute stress tolerance indices
#'
#' Compute several stress tolerance indices.
#' @param traits The list of traits.
#' @param geno The name of the column that identifies the genotypes.
#' @param normal The identification for the normal group.
#' @param stress The identification for the stress group.
#' @param dfr The name of the data frame.
#' @details The names of the columns for the traits must follow the convention
#' \code{trait.normal} and \code{trait.stress}, so the data frame must be in
#' wide format for the factor with levels \code{normal} and \code{stress}.
#' For a trait \code{y} with values \code{yn} and \code{ys} under normal and
#' stress conditions, the computed indices are:
#' \itemize{
#'  \item \code{tol} : Tolerance.
#'  \item \code{yrr} : Yield reduction ratio.
#'  \item \code{ssi} : Stress susceptibility index.
#'  \item \code{sti} : Stress tolerance index.
#'  \item \code{mp}  : Mean productivity.
#'  \item \code{gmp} : Geometric mean productivity.
#'  }
#' @return It returns a data frame with the indices and a data frame with
#' the stress intensity values for each trait. The names for the indices
#' follow the convention \code{trait.index}.
#' @author Raul Eyzaguirre.
#' @examples
#' traits <- c("nmtp", "mtwp", "nnomtp")
#' sti(traits, 'genotype', 'DTWW', 'DTWS', potatostress)
#' @export

sti <- function(traits, geno, normal, stress, dfr) {
  
  # Trait names
  
  traits.normal <- paste(traits, normal, sep = '.')
  traits.stress <- paste(traits, stress, sep = '.')
  
  # Stress intensity
  
  si.out <- NULL
  si.name <- NULL
  
  # Indices 
  
  for (i in 1:length(traits)) {
    
    # Normal value of the trait

    nv <- dfr[, traits.normal[i]]
    
    # Stress value of the trait
    
    sv <- dfr[, traits.stress[i]]
    
    # Tolerance
    
    index.name <- paste(traits[i], 'tol', sep = '.')
    dfr[, index.name] <- nv - sv
    
    # Yield reduction ratio
    
    index.name <- paste(traits[i], 'yrr', sep = '.')
    dfr[, index.name] <- 1 - sv / nv
    tmp <- index.name # Keep name for next index
    
    # Stress intensity
    
    si <- 1 - mean(sv) / mean(nv)
    si.out <- c(si.out, si)
    si.name <- c(si.name, traits[i])
    
    # Stress susceptibility index
    
    index.name <- paste(traits[i], 'ssi', sep = '.')
    dfr[, index.name] <- dfr[, tmp] / si
    
    # Stress tolerance index
    
    index.name <- paste(traits[i], 'sti', sep = '.')
    dfr[, index.name] <- nv * sv / mean(nv)^2

    # Mean productivity
    
    index.name <- paste(traits[i], 'mp', sep = '.')
    dfr[, index.name] <- (nv + sv) / 2

    # Geometric mean productivity
    
    index.name <- paste(traits[i], 'gmp', sep = '.')
    dfr[, index.name] <- sqrt(nv * sv)
   
  }
  
  # Return
  
  si.values <- data.frame(trait = si.name,
                          si = si.out)
    
  list(index.dfr = dfr,
       si.values = si.values)
  
}
