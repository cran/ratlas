% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatting.R
\name{formatting}
\alias{formatting}
\alias{fmt_count}
\alias{fmt_digits}
\alias{fmt_leading_zero}
\alias{fmt_minus}
\alias{fmt_replace_na}
\alias{fmt_corr}
\alias{fmt_prop}
\alias{fmt_prop_pct}
\title{Text and Number Formatting}
\usage{
fmt_count(x, big_interval = 3L, big_mark = ",")

fmt_digits(
  x,
  digits = 3,
  fmt_small = FALSE,
  max_value = NULL,
  keep_zero = FALSE
)

fmt_leading_zero(x)

fmt_minus(x, output = NULL)

fmt_replace_na(x, replacement = "&mdash;")

fmt_corr(x, digits, output = NULL)

fmt_prop(x, digits, fmt_small = TRUE, keep_zero = FALSE)

fmt_prop_pct(x, digits = 0, fmt_small = TRUE)
}
\arguments{
\item{x}{Number or number string to be formatted}

\item{big_interval}{Interval indicating where to place numeric dividers}

\item{big_mark}{Character used as mark between big interval before the
decimal}

\item{digits}{Number of decimal places to retain}

\item{fmt_small}{Indicator for replacing zero with \code{<} (e.g., \code{.000} becomes
\verb{<.001}). Default is \code{TRUE}.}

\item{max_value}{If \code{fmt_small} is \code{TRUE} and a \verb{max_value is supplied},
any value greater than the \code{max_value} is replaced with \code{>}
(e.g., if \code{max_value} = 50, then \code{60} becomes \verb{>49.9}). The number of digits
depends on \code{digits}.}

\item{keep_zero}{If \code{fmt_small} is \code{TRUE}, whether to preserve true 0s (e.g.,
\code{0.0000001} becomes \verb{<.001}, but \code{0.0000000} stays \code{.000}).}

\item{output}{The output type for the rendered document. One of \code{"latex"} or
\code{"html"}.}

\item{replacement}{The value to use when replacing missing values}
}
\value{
The updated character object of the same size as \code{x}.
}
\description{
These formatting functions are used to format numerical values in a
consistent manner. This is useful for printing numbers inline with text, as
well as for formatting tables. Many of the included formatting functions were
adapted from TJ Mahr's \href{https://github.com/tjmahr/printy}{printy} package.
}
\details{
\code{fmt_count()} is a wrapper for \code{\link[base:formatc]{base::prettyNum()}}. Prints a number with
a \code{big_mark} between every \code{big_interval}.

\code{fmt_digits()} is a wrapper for \code{\link[base:sprintf]{base::sprintf()}}. Prints a number with
\code{digits} number of decimal places, without losing trailing zeros, as happens
with \code{\link[base:Round]{base::round()}}.

\code{fmt_leading_zero()} removes the leading zero for decimal values.

\code{fmt_minus()} replaces hyphens with the HTML minus sign (\verb{&minus;}).

\code{fmt_replace_na()} replaces \code{NA} values with a specified replacement. This is
useful for formatting tables, when blanks are not desired. The default
behavior is to replace missing values with an em-dash (\verb{&mdash;}).

\code{fmt_prop_pct()} formats proportions as percentages. This takes a number
bounded between 0 and 1, multiplies it by 100, and then rounds to the
specified number of digits using \code{fmt_digits()}.

Two additional formatters are provided to format numbers according to the
American Psychological Association (APA) style guide. The 7th edition of the
APA style guide specifies that numbers bounded between [-1, 1] should not
include the leading zero (section 6.36; APA, 2020). This is the case for many
types of numbers commonly used by ATLAS including correlations, proportions,
probabilities, and \emph{p}-values. The \code{fmt_corr()} function is used to format
values bounded between [-1, 1]. Digits are first rounded to the specified
number of digits using \code{fmt_digits()}, and then leading zeros are removed
using \code{fmt_leading_zero()} and negative signs are replaced with
\code{fmt_minus()}. The \code{fmt_prop} is very similar, but is intended for values
between [0, 1]. This function also wraps \code{fmt_digits()} and
\code{fmt_leading_zero()}. However, \code{fmt_prop()} also replaces small values to
avoid values of \code{0} (e.g., \code{.00} is replaced with \verb{< .01}).
}
\examples{
test_cor <- cor(mtcars[, 1:4])
as.character(round(test_cor[1:4, 3], 2))
fmt_digits(test_cor[1:4, 3], 2)

fmt_digits(test_cor[1:4, 3], 2) \%>\%
  fmt_leading_zero()

fmt_digits(test_cor[1:4, 3], 2) \%>\%
  fmt_minus()

fmt_digits(c(test_cor[1:4, 3], NA_real_), 2) \%>\%
  fmt_replace_na(replacement = "&mdash;")

fmt_corr(test_cor[1:4, 3], 2)

fmt_prop(c(0.001, 0.035, 0.683), digits = 2)

}
\references{
American Psychological Association. (2020). \emph{Publication manual
of the American Psychological Association} (7th ed.).
\doi{doi:10.1037/0000165-000}
}
\seealso{
Other formatters: 
\code{\link{fmt_table}()},
\code{\link{padding}}
}
\concept{formatters}
