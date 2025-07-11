% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set-theme.R
\name{set_theme}
\alias{set_theme}
\title{Set default ggplot2 theme}
\usage{
set_theme(
  font = "Arial Narrow",
  discrete = c("okabeito", "atlas", "ggplot2"),
  continuous = c("viridis", "magma", "inferno", "plasma", "cividis", "ggplot2"),
  ...
)
}
\arguments{
\item{font}{The base font family to be used in plots.}

\item{discrete}{Color palette for discrete colors. One of "okabeito"
(default), "atlas", or "ggplot2".}

\item{continuous}{Color palette for continuous scales. One of "magma",
"inferno", "plasma", "viridis" (default), or "cividis", or "ggplot2".}

\item{...}{Additional arguments to pass to theme functions.}
}
\value{
None. Called for side effects.
}
\description{
Sets the default color schemes, fonts, and theme for ggplot2 plots. The
default color scheme for continuous variables is the
\href{https://CRAN.R-project.org/package=viridis}{viridis}
color palette, and the default color scheme for discrete variables is the
\href{https://jfly.uni-koeln.de/color/}{Okabe Ito} palette.
}
\examples{
\donttest{
set_theme("Arial Narrow")
}
}
