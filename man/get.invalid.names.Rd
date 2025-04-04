% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_names.R
\name{get.invalid.names}
\alias{get.invalid.names}
\title{Get column names not defined in crop ontology}
\usage{
get.invalid.names(dfr, add = NULL, crop = c("auto", "pt", "sp"))
}
\arguments{
\item{dfr}{The name of the data frame.}

\item{add}{Additional variables.}

\item{crop}{\code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.}
}
\value{
A character vector of invalid column names
}
\description{
Run \code{get_invalid_names()} after running \code{check.names()}
}
\details{
Check that fieldbook factors and variables' names correspond with the names defined
in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
CIP protocols.

Type \code{ptont()} or \code{spont()} to see the list of variables and
corresponding short labels and CO numbers.
}
\examples{
\dontrun{
tmp <- check.names(potatoyield)
get.invalid.names(tmp)
tmp <- check.names(pjpz09)
get.invalid.names(tmp)
}
}
\author{
Raul Eyzaguirre.
}
