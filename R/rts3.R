#' Response to selection with several locations and years
#'
#' It finds the optimum number of replications to get the maximum response
#' to selection with several locations and years for a given plot capacity,
#' number of locations, number of years, number of selected genotypes, genotypic variance,
#' genotypic by location variance, genotypic by year variance, genotypic by location by
#' year variance, and error variance.
#' @author Raul Eyzaguirre.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts3()} in the R console to run the app.
#' @return It returns a plot of response to seleccion versus number of replications
#' and computes the optimum number of replications and the response to selection
#' at this optimum value.
#' @export

rts3 <- function () {
  dirfiles <- paste(system.file(package = "st4gi"), "/shinyapps/rts3", sep="")
  shiny::runApp(dirfiles)
}
