% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_data_pt.R
\name{check.data.pt}
\alias{check.data.pt}
\title{Check consistency for potato experimental data}
\usage{
check.data.pt(
  dfr,
  f = 5,
  out.mod = c("none", "rcbd", "met"),
  out.max = 4,
  add = NULL,
  print.text = TRUE
)
}
\arguments{
\item{dfr}{The name of the data frame.}

\item{f}{Factor for extreme values detection. See details.}

\item{out.mod}{Statistical model for outliers' detection. See details.}

\item{out.max}{Threshold for outliers' detection.}

\item{add}{Additional quantitative traits.}

\item{print.text}{Logical, if \code{TRUE} the output is printed on screen.}
}
\value{
It returns:
\itemize{
\item \code{$Inconsist.List}, a \code{data.frame} with a list of all the
rows with some kind of inconsistency.
\item \code{$Inconsist.Matrix}, a \code{data.frame} with the positions
in the fieldbook data frame where inconsistencies occur. These are coded
with: (1) for inconsistencies among traits, (2) for out of range values,
(3) for extreme values or outliers.
}
}
\description{
Set of rules to check for consistency of potato experimental data.
Data labels must be defined as specified in
\url{http://www.cropontology.org/ontology/CO_330/Potato}.
}
\details{
The data frame must use the labels (lower or upper case) listed in
function \code{check.names.pt}.

Extreme low and high values are detected using the interquartile range.
The rule is to detect any value out of the interval 
\eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 5}.
If \code{f = 0}, the detection of extreme values is not executed.

Outliers are detected based on standardized residuals for some statistical
models. Options are \code{"rcbd"} and \code{"met"} for a randomized complete
block design and a multi environment trial with RCBD in each environment.
By default the threshold value is \code{out.max = 4}.
}
\author{
Raul Eyzaguirre.
}
