% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{inc}
\alias{inc}
\title{Generate a section for the yaml input}
\usage{
inc(input, sep = "\\n\\n  ")
}
\arguments{
\item{input}{a file containing markdown text}

\item{sep}{a separator for each line.}
}
\value{
a string
}
\description{
Generate a section for the yaml input
}
\examples{
\dontrun{
inc("front-matter/preface.Rmd")
}
}
