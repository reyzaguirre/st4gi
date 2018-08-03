#' Response to selection for a single experiment
#'
#' It finds the optimum number of replications to get the maximum response
#' to selection for a single experiment for a given plot capacity,
#' number of selected genotypes, genotypic variance and error variance.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts1()} in the R console to run the app.
#' @return It returns a plot of response to seleccion versus number of replications
#' and computes the optimum number of replications and the response to selection
#' at this optimum value.
#' @author Raul Eyzaguirre.
#' @export

rts1 <- function () {
  
  dirfiles <- paste0(system.file(package = "st4gi"), "/shinyapps/rts1")
  shiny::runApp(dirfiles)

}

#' Response to selection with several locations
#'
#' It finds the optimum number of replications to get the maximum response
#' to selection with several locations for a given plot capacity,
#' number of locations, number of selected genotypes, genotypic variance,
#' genotypic by location variance, and error variance.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts2()} in the R console to run the app.
#' @return It returns a plot of response to seleccion versus number of replications
#' and computes the optimum number of replications and the response to selection
#' at this optimum value.
#' @author Raul Eyzaguirre.
#' @export

rts2 <- function () {
  
  dirfiles <- paste0(system.file(package = "st4gi"), "/shinyapps/rts2")
  shiny::runApp(dirfiles)
  
}

#' Response to selection with several locations and years
#'
#' It finds the optimum number of replications to get the maximum response
#' to selection with several locations and years for a given plot capacity,
#' number of locations, number of years, number of selected genotypes, genotypic variance,
#' genotypic by location variance, genotypic by year variance, genotypic by location by
#' year variance, and error variance.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts3()} in the R console to run the app.
#' @return It returns a plot of response to seleccion versus number of replications
#' and computes the optimum number of replications and the response to selection
#' at this optimum value.
#' @author Raul Eyzaguirre.
#' @export

rts3 <- function () {
  
  dirfiles <- paste0(system.file(package = "st4gi"), "/shinyapps/rts3")
  shiny::runApp(dirfiles)
  
}

#' Response to selection with several locations in two steps (two years)
#'
#' It computes the response to selection for each step in a two steps selection
#' with several locations for a given number of genotypes at step 1,
#' the number of locations, replications and selected genotypes at step 1,
#' the number of locations, replications and selected genotypes at step 2,
#' and the genotypic, genotypic by location, genotypic by year,
#' genotypic by location by year, and error variances.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts4()} in the R console to run the app.
#' @return It returns the response to selection at step 1 and 2.
#' @author Raul Eyzaguirre.
#' @export

rts4 <- function () {
  
  dirfiles <- paste0(system.file(package = "st4gi"), "/shinyapps/rts4")
  shiny::runApp(dirfiles)
  
}
