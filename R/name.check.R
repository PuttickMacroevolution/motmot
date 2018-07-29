#' Name check (internal function)
#'
#' This is an internal function to check if names in a trait vector or matrix match a tree
#' @param phy a phylogeny in \pkg{ape} phylo format
#' @param data A matrix or vector of named trait data corresponding to species in phy
#' @author Gavin Thomas
#' @export

name.check <- function (phy, data) 
{
    if (is.vector(data)) {
        data.names <- names(data)
    }
    else {
        data.names <- rownames(data)
    }
    t <- phy$tip.label
    r1 <- t[is.na(match(t, data.names))]
    r2 <- data.names[is.na(match(data.names, t))]
    r <- list(sort(r1), sort(r2))
    names(r) <- cbind("tree_not_data", "data_not_tree")
    if (length(r1) == 0 && length(r2) == 0) 
	return("OK")
    else return(r)
}