% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_pos.R
\name{check.pos}
\alias{check.pos}
\title{Check row and column positions}
\usage{
check.pos(dfr, row = "row", col = "col", rep = NULL)
}
\arguments{
\item{dfr}{The name of the data frame.}

\item{row}{The name of the column that identifies the rows.}

\item{col}{The name of the column that identifies the columns.}

\item{rep}{The name of the column that identifies the replications.}
}
\value{
A list of plots (unique row and column position)
with more than one genotype.
}
\description{
This function checks that there is only one genotype in each row and column position.
}
\details{
If \code{rep} is not \code{NULL}, then it checks
positions for each replication.
}
\examples{
# Create a design
dfr <- cr.rcbd(1:20, 3, 10)
dfr <- dfr$book
# Check positions
check.pos(dfr, 'row', 'col', 'block')
}
\author{
Raul Eyzaguirre.
}
