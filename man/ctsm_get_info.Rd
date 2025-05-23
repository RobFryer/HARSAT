% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information_functions.R
\name{ctsm_get_info}
\alias{ctsm_get_info}
\title{Gets data from the standard reference tables}
\usage{
ctsm_get_info(
  ref_table,
  input,
  output,
  compartment = NULL,
  na_action = c("fail", "input_ok", "output_ok", "ok"),
  info_type = NULL,
  sep = "."
)
}
\arguments{
\item{ref_table}{the reference table}

\item{input}{the input}

\item{output}{the output}

\item{compartment}{the compartment}

\item{na_action}{character, how to handle missing values, one of:
\code{fail} (no missing values allowed), \code{input_ok} (allows missing values in input, but all
non-missing values must be recognised, and must have output),
\verb{output_ok`` requires all input values to be present, but allows missing values in  output (for e.g. dryweight by species), or }ok` (allows missing values everywhere)}

\item{info_type}{the information type, by default picked up from the function call}

\item{sep}{character, the separator, defaulting to '.'}
}
\description{
Gets data from the standard reference tables
}
