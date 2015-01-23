#'  Elston Index
#'
#' Function to compute the Elston index (Elston, R. C., 1963).
#' @param traits List of traits.
#' @param geno The genotypes.
#' @param data The name of the data frame containing the data.
#' @param lb Lower bound. 1 for k = min(x) and 2 for k = (n*min(x) - max(x))/(n-1)
#' @author Raul Eyzaguirre
#' @details The Elston index is a weight free index.
#' @return It returns the Elston index value and the Elston index value
#' sorted in descending order.
#' @references
#' Elston, R. C. (1963). A weight-free index for the purpose of ranking or selection
#' with respect to several traits at a time. Biometrics. 19(1): 85-97.
#' @examples
#' # The data
#' head(spg)
#' str(spg)
#'
#' # Run Elston index with all the traits
#' elston(c("rytha", "bc", "dm", "star", "nocr"), "geno", spg)
#' @export

elston <- function(traits, geno, data, lb = 1) {

  # inits

  nt <- length(traits) # number of traits
  k <- NULL
  ng <- nlevels(factor(data[,geno])) # number of genotypes

  # compute standardized means

  df <- data.frame(tapply(data[,traits[1]], data[,geno], mean, na.rm = T))
  if (nt > 1){
    for (i in 2:nt)
      df <- cbind(df, tapply(data[,traits[i]], data[,geno], mean, na.rm = T))
    for (i in 1:nt)
      df[,i+nt] <- (df[,i] - mean(df[,i]))/sd(df[,i])
  }

  # compute lower bounds

  if (lb == 1)
    for (i in 1:nt)
      k[i] <- min(df[,nt+i])

  if (lb == 2)
    for (i in 1:nt)
      k[i] <- (ng * min(df[,nt+i]) - max(df[,nt+i]))/(ng-1)

  # Elston index

  df$EI <- df[,nt+1] - k[1]
  if (nt > 1)
    for (i in 2:nt)
      df$EI <- df$EI * (df[,nt+i] - k[i])

  orden <- order(df$EI, decreasing = T)

  # results

  list(Elston.Index = df$EI, Sorted.Elston.Index = df$EI[orden])
}
