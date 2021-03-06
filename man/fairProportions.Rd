% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fairProportions.R
\name{fairProportions}
\alias{fairProportions}
\title{Calculate fair proportions phylogenetic diversity metric}
\usage{
fairProportions(phy, nodeCount = FALSE)
}
\arguments{
\item{phy}{An object of class \code{phylo} (see \pkg{ape}).}

\item{nodeCount}{Logical - should root to tip node counts be returned (default is \code{FALSE})}
}
\value{
Returns a matrix of fair proportion for all tips in phylogeny and node counts if selected.
}
\description{
Calculate fair proportions phylogenetic diversity metric
Note that \code{as.rateMatrix} calls the CAIC function \code{vcv.array} multiple times and this can be slow for large phylogenies (though faster than using the "ape" equivalent \code{vcv.phylo}).
}
\examples{
data(anolis.tree)

fp <- fairProportions(anolis.tree)
fpNodes <- fairProportions(anolis.tree, nodeCount=TRUE)
}
\references{
Redding, D.W. and Mooers, A.O. (2006). Incorporating evolutionary measures into conservation prioritisation. Conservation Biology, 20, 1670-1678.

Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
}
\author{
Gavin Thomas
}
