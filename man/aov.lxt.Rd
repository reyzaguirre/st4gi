% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aov_lxt.R
\name{aov.lxt}
\alias{aov.lxt}
\title{Analysis for a Line x Tester genetic design.}
\usage{
aov.lxt(dfr, y, line, tester, rep)
}
\arguments{
\item{dfr}{The name of the data frame.}

\item{y}{The name of the column for the variable to analyze.}

\item{line}{The name of the column that identifies the lines.}

\item{tester}{The name of the column that identifies the testers.}

\item{rep}{The name of the column that identifies the replications.}
}
\value{
It returns the analysis for the line by tester genetic design.
}
\description{
This is a modified version of the \code{"lineXtester"} function of package
\code{"agricolae"} for the analysis of a Lina x Tester genetic design.
}
\details{
Data must correspond with a line by tester design and must be balanced.
If there are missing values they are not estimated and results can be dubious.
}
\examples{
aov.lxt(lxt, "yield", "line", "tester", "rep")
}
\author{
Raul Eyzaguirre.
}
