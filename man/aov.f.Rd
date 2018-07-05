% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aovf.R
\name{aov.f}
\alias{aov.f}
\title{ANOVA for a factorial experiment}
\usage{
aov.f(trait, factors, rep, design = c("crd", "rcbd"), data, maxp = 0.1)
}
\arguments{
\item{trait}{The trait to analyze.}

\item{factors}{The factors.}

\item{rep}{The replications or blocks.}

\item{design}{The statistical design, \code{crd} or \code{rcbd}.}

\item{data}{The name of the data frame.}

\item{maxp}{Maximum allowed proportion of missing values to estimate, default is 10\%.}
}
\value{
It returns the ANOVA table.
}
\description{
Fit an analysis of variance model for a factorial experiment with a CRD or RCBD
}
\details{
If data is unbalanced, missing values are estimated up to an specified maximum
proportion, 10\% by default.
}
\examples{
aov.f("asc.dw", c("geno", "treat"), "rep", "crd", asc)
}
\author{
Raul Eyzaguirre.
}