% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render.R
\name{techreport_pdf}
\alias{techreport_pdf}
\title{Create an R Markdown PDF Document Tech Report}
\usage{
techreport_pdf(apa6 = FALSE, ...)
}
\arguments{
\item{apa6}{Should the old}

\item{...}{Arguments to be passed to \verb{[bookdown::pdf_document2]}}
}
\value{
A modified \code{pdf_document2} with the standard tech report formatting.
}
\description{
This is a function called in the output of the yaml of the Rmd file to
specify using the standard DLM tech report pdf document formatting.
}
\examples{
\dontrun{
output: ratlas::techreport_pdf
}
}
