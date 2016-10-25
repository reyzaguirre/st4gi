#' Data transformations
#'
#' This function performs different data transformations.
#' @param trait The trait to transform.
#' @param type The transformation type. See details.
#' @param base Base for the logarithmic transformation. Base 10 by default.
#' @param n Additional parameter for arc-sine transformation. See details.
#' @param data The name of the data frame containing the data.
#' @author Raul Eyzaguirre.
#' @details Available transformations are:
#' 
#' \code{none} for no transformation.
#' 
#' \code{logy} for the logaritmic transformation log(y). This transformation is
#' recomended for data that follow a multiplicative intead of an additive model.
#' 
#' \code{logy1} for the logaritmic transformation log(y + 1). The same as the
#' previous case, but when the data set includes small values (e.g. less than 10).
#' 
#' \code{sqrty} for the square root transformation sqrt(y). This transformation is
#' recomended for count data, which typically follow a Poisson distribution where
#' the variance is proportional to the mean. It is also recommended for percentage
#' data where the range is between 0 and 20% or between 80 and 100%. However note
#' that for Poisson data, a Poisson regression model could be a better option.
#' 
#' \code{sqrty1} for the square root transformation sqrt(y + 0.5). The same as the
#' previous case but when most of the values in the data set are less than 10,
#' especially if zeros are present.
#' 
#' \code{arcsin} for the  arc-sine transformation arcsin(y^0.5). This transformation
#' is recomended for data on proportions, which typically follow a binomial
#' distribution. Where the values of 0 or 1 are present, these must be substituted
#' by 1/4n and 1-1/4n, where \code{n} is the denomitator for the computation of the
#' proportions. When the proportions are in the range 0.2 to 0.8 no transformation
#' could be needed, and where some are on either the range 0 to 0.2 or 0.8 to 1 a
#' square root transformation could be useful. Finally, Note that for binomial data,
#' a binomial regression model could be a better option.
#' @return It returns the transformed trait.
#' @examples
#' dtr("nonc", "logy", data = pjpz09)
#' @export

dtr <- function(trait, type = c("none", "logy", "logy1", "sqrty", "sqrty1", "arcsin"),
                base = 10, n = NULL, data) {

  # match arguments
  
  type <- match.arg(type)
  
  # log transformations
  
  if (type == "logy") {
    if (sum(data[, trait] <= 0) > 0) {
      data[data[, trait] <= 0, trait] <- NA
      warning("Values <= 0 converted to NA", call. = FALSE)
    }
    data[, trait] <- log(data[, trait], base)
  }

  if (type == "logy1"){
    if (sum(data[, trait] + 1 <= 0) > 0) {
      data[data[, trait] + 1 <= 0, trait] <- NA
      warning("Values <= -1 converted to NA", call. = FALSE)
    }
    data[, trait] <- log(data[, trait] + 1, base)
  }
  
  # sqrt transformation

  if (type == "sqrty") {
    if (sum(data[, trait] < 0) > 0) {
      data[data[, trait] < 0, trait] <- NA
      warning("Values < 0 converted to NA", call. = FALSE)
    }
    data[, trait] <- sqrt(data[, trait])
  }
  
  if (type == "sqrty1") {
    if (sum(data[, trait] + 0.5 < 0) > 0) {
      data[data[, trait] + 0.5 < 0, trait] <- NA
      warning("Values < -0.5 converted to NA", call. = FALSE)
    }
    data[, trait] <- sqrt(data[, trait] + 0.5)
  }
  
  # arc-sine transformation
  
  if (type == "arcsin") {
    if (sum(data[, trait] < 0) > 0 | sum(data[, trait] > 1) > 0) {
      data[data[, trait] < 0, trait] <- NA
      data[data[, trait] > 1, trait] <- NA
      warning("Values < 0 or > 1 converted to NA", call. = FALSE)
    }
    if (sum(data[, trait] == 0) > 0 | sum(data[, trait] == 1) > 0) {
      if (!is.null(n)) {
        data[data[, trait] == 0, trait] <- 1/(4*n)
        data[data[, trait] == 1, trait] <- 1 - 1/(4*n)
        warning("Values = 0 and = 1 replaced with 1/4n and 1 - 1/4n", call. = FALSE)
      }
      if (is.null(n)) {
        warning("n is not specified. Values = 0 or = 1 are not replaced", call. = FALSE)
      }
    }
    data[, trait] <- sqrt(data[, trait] + 0.5)
  }

  # results
  
  data
}