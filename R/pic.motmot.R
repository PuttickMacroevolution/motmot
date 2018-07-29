#' @title Phylogenetically independent contrasts (internal)
#' @description Calculates phylogenetically independent contrasts.
#' @param x A matrix of trait values with taxon names as rownames.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @details Extracts values for contrasts, expected variances using contrasts by calling \code{pic}.
#' @return contr A matrix with two columns containing raw contrasts in the first column and their expected variances in the second column.
#' @return root.v Expected variances of branches either side of root
#' @return V Expected variance at the root 
#' @references Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
#' @author Gavin Thomas, Rob Freckleton, Emmanuel Paradis
#' @export

pic.motmot <- function (x, phy) {
out <- pic(x, phy, rescaled.tree=TRUE, var.contrasts=TRUE, scaled=FALSE)
   contr <- out[[1]]	
	nb.tip <- length(phy$tip.label)
	idx <- which(out[[2]]$edge[,1] == (nb.tip+1) )
	
	root.v <- out[[2]]$edge.length[idx]
	V <- prod(root.v)/(sum(root.v))
  
  return(list(contr = contr, root.v = root.v, V = V))
}
