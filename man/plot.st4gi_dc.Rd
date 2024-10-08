% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_checkdata.R
\name{plot.st4gi_dc}
\alias{plot.st4gi_dc}
\title{Data consistency visualization}
\usage{
\method{plot}{st4gi_dc}(x, ...)
}
\arguments{
\item{x}{An object of class \code{checkdata}.}

\item{...}{Additional plot arguments.}
}
\value{
It returns a plot.
}
\description{
This function produces a visualization of data consistency.
}
\details{
It produces a visual matrix representation of data consistency problems.
}
\examples{
checks <- check.data(pjpz09)
plot(checks)
checks <- check.data(potatoyield)
plot(checks)
}
\author{
Raul Eyzaguirre.
}
