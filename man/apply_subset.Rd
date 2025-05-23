% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics_functions.R
\name{apply_subset}
\alias{apply_subset}
\title{Apply subsetting to a time series}
\usage{
apply_subset(timeSeries, subset, env = parent.frame())
}
\arguments{
\item{timeSeries}{a time series data frame}

\item{subset}{(default \code{NULL}) either \code{NULL}, which selects all entries, or a logical,
which \code{TRUE} selects all entries and \code{FALSE} none of them, or an expression, which
will be evaluated in the context of the dataframe to generate a subsetting vector.}

\item{env}{the calling environment for variable values.}
}
\description{
This is an internal function that applies a subsetting function to a timeseries.
It is somewhat complex, due to way it is designed to be called. The subset is designed
to be passed either as a value (\code{NULL} or  logical), as a expression, or as a variable
holding an expression. When it is an expression, it is applied in the context of the
timeseries data frame to generate a vector of booleans for subsetting. The complexity
comes from the way this implements lazy evaluation, so we cannot evaluate the expression
in the normal calling context.

The function also removes row names, and converts the series column to row names.
}
