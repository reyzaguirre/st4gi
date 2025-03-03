% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ck_abd.R
\name{ck.abd}
\alias{ck.abd}
\title{Check data for an ABD}
\usage{
ck.abd(dfr, y, geno, rep)
}
\arguments{
\item{dfr}{The name of the data frame.}

\item{y}{The name of the column for the variable to analyze.}

\item{geno}{The name of the column that identifies the genotypes including checks.}

\item{rep}{The name of the column that identifies the replications.}
}
\value{
The number of checks \code{ng.ck}, the number of no checks \code{ng},
the number of missing values for checks \code{nmis.ck}), the number of missing
values for no checks \code{nmis}, the number \code{nck.0} and list \code{ck.0}
of checks without data, the number \code{nck.1} and list \code{ck.1} of checks
with only one datum, the number of checks with at least two data \code{nck.2},
the number of checks that appear more than once in a given replication
(\code{nck.mult}), the number of replications \code{nrep}, and the number of rows
in the data frame with missing values for factors (\code{nmis.fac}).
}
\description{
This function checks the frequencies of genotypes in an ABD.
}
\examples{
# Create design
dfr <- cr.abd(1:50, c('a', 'b', 'd'), 5, 10)
dfr <- dfr$book
# Create some random data
dfr$y <- rnorm(65)
# Delete some values
dfr[c(1, 5, 7, 56), 'y'] <- NA
# Delete some values for classification factors
dfr[64, 'geno'] <- NA
# Check the design
ck.abd(dfr, 'y', 'geno', 'block')
}
\author{
Raul Eyzaguirre.
}
