#' @title Drop tips from a phylogenetic tree while preserving deleted nodes
#' @description Wrapper for the \code{ape} function \code{drop.tip} that preserves the number of nodes affecting each branch. For use with the \code{psi} and \code{multipsi} models.
#' @param phy Phylogenetic tree in \code{phylo} format 
#' @param tip A vector of mode numeric or character specifying the tips to delete, to be passed to \code{drop.tip}
#' @return Phylogenetic tree in \code{phylo} format, with an added element \code{Shid}, a vector of numbers of observed but "missing" speciation events per branch, in the same order as the branches in the \code{phylo} object
#' @references Ingram, T. 2011. Speciation along a depth gradient in a marine adaptive radiation. Proc. R. Soc. B 278: 613-618.
#' @author Travis Ingram
#' @seealso \code{\link{transformPhylo.ML}}
#' @export

dropTipSobs<-function(phy, tip)
{
	phy1 <- phy
	phy1$edge.length <- phy$edge.length^0
	phy1 <- drop.tip(phy1, tip)
	phy2 <- drop.tip(phy, tip)
	phy2$Sobs <- phy1$edge.length
	phy2
}
