#' allCladeMembers (internal function)
#'
#' This is an internal function to generate an \code{allCladeMembersMatrix}
#' @param phy phylogeny in \pkg{ape} \code{phylo} format
#' @return cladeMembersMatrix
#' @author Gavin Thomas
#' @export

allCladeMembers <- function(phy) {
	
	nodeIDs <- c(1:Ntip(phy), (Ntip(phy)+2):(Ntip(phy) + Nnode(phy)))
	k <- length(nodeIDs)
   cladeMembersMatrix <-  sapply(nodeIDs, function(k) {
    	nodeShiftID <- c(k, node.descendents(x = k , phy = phy))
    	as.numeric(phy$edge[, 2] %in% nodeShiftID)
    	}
    )
    return(cladeMembersMatrix)
}