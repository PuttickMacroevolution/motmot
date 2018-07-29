#' @title Sample hidden speciation events along branches of a tree (internal function)
#' @description Uses estimated speciation and extinction rates to sample the number of speciation events 'hidden' by subsequent extinction on each branch of a tree following Bokma (2008). For use with the \code{psi} and \code{multipsi} models. 
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param lambda.sp Estimate of the rate of speciation "lambda"
#' @param mu.ext Estimate of the rate of extinction "mu"
#' @param useMean A logical indicating whether to output the average or expected number of hidden speciation events per branch, which may be non-integer (if \code{TRUE}), or to sample an integer number on each branch from a Poisson distribution (if \code{FALSE}, the default)
#' @details The expected number of hidden speciation events are calculated for each branch given its start and end times, and estimates of lambda and mu which are assumed to be constant across the tree. To properly account for uncertainty in the effect of extinction on the number of nodes affecting each branch of a tree, it may be appropriate to repeat model-fitting on many realizations of \code{Sobs} on the tree of interest (similar to evaluating phylogenetic uncertainty)
#' @return Phylogenetic tree in \code{phylo} format, with an added element \code{Sobs}, a vector of numbers of hidden speciation events per branch, in the same order as the branches in the \code{phylo} object
#' @seealso \code{\link{transformPhylo.ML}}
#' @references Bokma, F. 2008. Detection of "punctuated equilibrium" by Bayesian estimation of speciation and extinction rates, ancestral character states, and rates of anagenetic and cladogenetic evolution on a molecular phylogeny. Evolution 62: 2718-2726.
#' @references Ingram, T. 2011. Speciation along a depth gradient in a marine adaptive radiation. Proc. R. Soc. B 278: 613-618.
#' @author Travis Ingram
#' @export

sampleHiddenSp <- function (phy, lambda.sp = NULL, mu.ext = NULL, useMean = FALSE) 
{
    if (is.null(lambda.sp) | is.null(mu.ext)) stop("Please provide values for lambda and mu ")  
    phy$node.label <- NULL
    edge.node.times <- nodeTimes(phy)
    to <- edge.node.times[ , 1]
    te <- edge.node.times[ , 2]
    te[te < 0] <- 0
    if (mu.ext > 0) {
        if (lambda.sp == mu.ext) {
            expSh <- 2 * lambda.sp * (to - te) + 2 * log((1 + lambda.sp * te) / (1 + lambda.sp * to))
        } else {
            expSh <- 2 * lambda.sp * (to - te) + 2 * log((lambda.sp * exp(te *  (lambda.sp - mu.ext)) - mu.ext) / (lambda.sp * exp(to * (lambda.sp - mu.ext)) - mu.ext))
        }
        if (useMean) {
            hidden.speciation <- expSh
        } else {
            hidden.speciation <- rpois(length(expSh), expSh)
        }
	} else {
        hidden.speciation <- rep(0, length(to))
    }
    phy$hidden.speciation <- hidden.speciation
    phy
}
