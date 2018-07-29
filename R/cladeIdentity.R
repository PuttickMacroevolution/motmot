#' @title Identify branches (including tips) descended from a node (internal function).
#' @description Internal function to get presence absence of descendent branches from a vector of node numbers. The descendents include the branch leading to the focal node (i.e. node defines the stem group no crown group
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param nodeIDs Vector of node numbers (positive integers).
#' @param cladeMembersObj Matrix of clade membership
#' @details The function returns a matrix of unique presences given the selected node. If the selected nodes are nested then presences are only recorded for the least inclusive node.
#' @return matrix Matrix of unique presences for each node id
#' @author Gavin Thomas
#' @examples
#' ## Read in phylogeny and data from Thomas et al. (2009)
#' data(anolis.tree)
#' data(anolis.data)
#'
#' cladeIdentityMatrix <- cladeIdentity(phy=anolis.tree, nodeIDs=170)
#' @export

cladeIdentity <- function (phy, nodeIDs, cladeMembersObj=NULL) 
{
	
	k <- length(nodeIDs)
	
    if(is.null(cladeMembersObj)) {
		cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 
														   1]))
		for (i in 1:k) {
			nodeShiftID <- c(nodeIDs[i], node.descendents(x = nodeIDs[i], phy = phy))
			cladeMembers[, i] <- as.numeric(phy$edge[, 2] %in% nodeShiftID)
		}
	}
	
	if (is.null(cladeMembersObj)==FALSE) {
		allNodes <- c(1:Ntip(phy), (Ntip(phy)+2):(length(phy$edge.length)+1))   #######
		cladeMembers <- as.matrix(cladeMembersObj[,match(nodeIDs, allNodes)])
	}
	
    originalOrder <- colSums(cladeMembers)
    richnessOrder <- sort(originalOrder, decreasing = FALSE, 
						  index.return = TRUE)
    cladeMembersOrdered <- matrix(cladeMembers[, richnessOrder$ix], 
								  ncol = length(nodeIDs))
	
	if (k>1) {
		for (i in 2:k){
			if (i ==2 ) {	cladeMembersOrdered[,i] <- cladeMembersOrdered[,i] - cladeMembersOrdered[,1:(i-1)]}
			else {cladeMembersOrdered[,i] <- cladeMembersOrdered[,i] - rowSums(cladeMembersOrdered[,1:(i-1)])}
		}
	}
	
    cladeMembers <- cladeMembersOrdered[, sort(richnessOrder$ix, index.return = TRUE)$ix]
    cladeMembers <- matrix(cladeMembers, ncol = length(nodeIDs))
    return(cladeMembers)
}

