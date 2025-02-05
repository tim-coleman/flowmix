% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-2d.R
\name{plot_2d_threepanels}
\alias{plot_2d_threepanels}
\title{Plot three 2d panels of data, optionally with a model (from \code{obj}).}
\usage{
plot_2d_threepanels(
  obj = NULL,
  ylist,
  countslist = NULL,
  tt,
  labels = NULL,
  plist_return = FALSE,
  colours = NULL,
  cruise_id = NULL
)
}
\arguments{
\item{obj}{A flowmix object.}

\item{ylist}{Data.}

\item{countslist}{Defaults to NULL.}

\item{tt}{time point of interest, out of 1 through \code{length(ylist)}.}
}
\value{
A grob object containing a 3-panel plot.
}
\description{
Plot three 2d panels of data, optionally with a model (from \code{obj}).
}
