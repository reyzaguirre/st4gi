#' Response to selection with several locations
#'
#' It finds the optimum number of replications to get the maximum response
#' to selection with several locations for a given plot capacity,
#' number of locations, number of selected genotypes, genotypic variance,
#' genotypic by location variance, and error variance.
#' @author Raul Eyzaguirre.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts2()} in the R console to run the app.
#' @return It returns a plot of response to seleccion versus number of replications
#' and computes the optimum number of replications and the response to selection
#' at this optimum value.
#' @export

rts2 <- function () {
  dirfiles <- paste(system.file(package = "st4gi"), "/shinyapps/rts2", sep="")
  shiny::runApp(dirfiles)
}
