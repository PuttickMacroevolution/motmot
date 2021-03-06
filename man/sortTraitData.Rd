% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sortTraitData.R
\name{sortTraitData}
\alias{sortTraitData}
\title{Sort data and remove missing entries for tree and trait data}
\usage{
sortTraitData(
  phy,
  y,
  data.name = NULL,
  log.trait = TRUE,
  pass.ultrametric = FALSE
)
}
\arguments{
\item{phy}{An object of class \code{phylo} or \code{multiPhylo}  (see \pkg{ape})}

\item{y}{A matrix of trait values with taxon names as rownames. Missing values should be NA}

\item{data.name}{If null the first column of y is assummed as the trait, otherwise if y is a matrix with more than one column either the name of the column or the number of the column must be supplied by data.name}

\item{log.trait}{Logical. If \code{TRUE}, data are log-transformed}

\item{pass.ultrametric}{Although trees that are believed to be ultrametric to pass the function \code{is.ultrametric} in \pkg{ape}}
}
\value{
phy Tree with missing data pruned

trait Rearranged data with missing species removed
}
\description{
Plots a phylogeny with lines representing the value of a continuous trait
}
\examples{
data(anolis.tree)
data(anolis.data)
attach(anolis.data)
male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
any(is.na(male.length[,1]))
data.sorted <- sortTraitData(anolis.tree, male.length)
phy <- data.sorted[[1]]
male.length <- data.sorted[[2]]
}
\author{
Mark Puttick
}
