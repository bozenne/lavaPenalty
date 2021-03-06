% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Penalty_simulation.R
\name{simForm}
\alias{simForm}
\title{Spatial simulation}
\usage{
simForm(n.obs, xmax, ymax, radius, center = NULL, coords.centered = TRUE,
  distance = "euclidean")
}
\arguments{
\item{n.obs}{the number of repetitions.}

\item{xmax}{the length of the spatial field.}

\item{ymax}{the width of the spatial field}

\item{radius}{the maximum distance to the center for which pixels will set to 1. The others will be set to 0.}

\item{center}{the position relative to which the distance will be computed.}

\item{coords.centered}{should the coordinates be centered around 0,0? Else they will start at 0,0 and take positve values. May affect the computation of the distance.}

\item{distance}{the type of distance to be used.}
}
\value{
a list containing:
\itemize{
\item{"coords"} the coordinates of each point on the field
\item{"center"} the center used to compute the distance
\item{"distCenter"} the distance of each point to the center
\item{"X"} the distance thresholded by the radius
}
}
\description{
Simulate a geometrical figure on a 2D spatial field
}
\details{
The available distances are those of the \code{dist} function from the \emph{stats} package.
}
\examples{
image(simForm(100, 10, 10, 2)$X)
image(simForm(100, 10, 10, 2, distance = "maximum")$X)

}
