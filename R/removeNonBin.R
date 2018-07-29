#' @title remove species occuring before time in the past (internal function) 
#'
#' @description removes tips and lineages after a time in the past
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param traitData data associated with the species
#' @param keepByTime an age at which to keep preserve tips before the present (time = 0)
#' @param tol edge length precision in time cut (default = 1e-08)
#' @return a list with the prunedPhylogeny and a prunedData
#' @references Puttick, M. N., Kriwet, J., Wen, W., Hu, S., Thomas, G. H., & Benton, M. J. (2017). Body length of bony fishes was not a selective factor during the biggest mass extinction of all time. Palaeontology, 60, 727-741.
#' @author Mark Puttick

removeNonBin <- function(phy, traitData, keepByTime=0, tol=1e-08){

    all.times <- nodeTimes(phy)
    tip.time <- all.times[which(phy$edge[,2] <= Ntip(phy)) , 2]
	non.bin <- which(tip.time > (keepByTime - tol))
	traitLess <- traitData[ - non.bin]
   	phylo <- drop.tip(phy, non.bin)
	return(list(prunedPhy=phylo, prunedData=traitLess))
}