#' Response to selection with several locations in two steps (two years)
#'
#' It computes the response to selection for each step in a two steps selection
#' with several locations for a given number of genotypes at step 1,
#' the number of locations, replications and selected genotypes at step 1,
#' the number of locations, replications and selected genotypes at step 2,
#' and the genotypic, genotypic by location, genotypic by year,
#' genotypic by location by year, and error variances.
#' @author Raul Eyzaguirre.
#' @details It uses package \code{shiny} for the web layout.
#' Type \code{rts4()} in the R console to run the app.
#' @return It returns the response to selection at step 1 and 2.
#' @export

rts4 <- function () {
  dirfiles <- paste(system.file(package = "st4gi"), "/shinyapps/rts4", sep="")
  shiny::runApp(dirfiles)
}
