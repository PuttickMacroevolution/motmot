#' Extimate ancestral state (internal function)
#'
#' This is an internal function to estimate the character value at the root node
#' @param phy a phylogeny in APE phylo format
#' @param y A matrix of trait data corresponding to species in phy
#' @author Gavin Thomas

ancState <- function (phy, y) {

			V <- vcv(phy)
			x <- rep(1, Ntip(phy))
			iV <- solve(V)
			xVix <- crossprod(x, iV %*% x)
			xViy <- crossprod(x, iV %*% y)
			mu <- solve(xVix) %*% xViy
			return(as.numeric(mu))
  }  
