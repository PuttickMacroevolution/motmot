% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloVar.R
\name{phyloVar}
\alias{phyloVar}
\title{Calculation of Brownian variance.}
\usage{
phyloVar(
  rateData,
  rate = NULL,
  common.mean = FALSE,
  lambda.est = TRUE,
  lambda = 1,
  meserr = FALSE
)
}
\arguments{
\item{rateData}{an object of class \code{rateData}}

\item{rate}{a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal.}

\item{common.mean}{a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.}

\item{lambda.est}{Logical. Fit Pagel's lambda.}

\item{lambda}{Numeric value for lambda from 0-1.}

\item{meserr}{Logical. Include measurement error.}
}
\value{
phylo.var  phylogenetic variance (Brownian variance)
}
\description{
This function calculates the phylogenetic variance (Brownian variance, or rate) of the data given the tree and model of evolution
}
\examples{
 ## Read in phylogeny and data from Thomas et al. (2009)
data(anolis.tree)
data(anolis.data)

## Convert data to class rateData with a rateMatrix object as input
anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)

anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)

# A model with a different rate in each of the four groups. The 'fixed' command is used to determine
# whether a particular rate is to be constrained or not. Use '1' to fix a group and 'FALSE' to show
# that the parameter is not fixed and should be estimated. The values should be entered in the same 
# order as the ranking of the groups. That is, group 0 (small islands) takes position one in the 
# fixed vector, group 1 (large island trunk crown and trunk ground) takes position 2 and so on. 
# The default is to allow each group to take a different mean. 

phyloVar(anolis.rateData, rate=c(1,2,1,1), common.mean=FALSE)
# common mean for all groups
phyloVar(anolis.rateData, rate=c(1,2,1,1), common.mean=TRUE)
}
\references{
Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
 @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
}
\author{
Gavin Thomas, Rob Freckleton
}
