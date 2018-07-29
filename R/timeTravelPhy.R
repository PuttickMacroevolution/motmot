#' @title timeTravelPhy (internal function) 
#'
#' @description removes tips and lineages after a time in the past
#' @param phy An object of class "phylo" (see ape package).
#' @param node nodes arising more recently than the cut time
#' @param nodeEstimate trait the number of descendants arising from the nodes
#' @param timeCut position at which to cut the phylogeny
#' @param traits Logical. Include trait values in the output tree
#' @return the pruned phylogeny and a 'tipObject' of the number of lineages found in the pruned branches
#' @references Puttick, M. N., Kriwet, J., Wen, W., Hu, S., Thomas, G. H., & Benton, M. J. (2017). Body length of bony fishes was not a selective factor during the biggest mass extinction of all time. Palaeontology, 60, 727-741.
#' @author Mark Puttick

timeTravelPhy <- function(phy, node, nodeEstimate, timeCut, traits=TRUE){

	countDesAll <- c()
	node.in <- match(unique(phy$edge[,1]), phy$edge[, 1])
	originalNodeTimes <- nodeTimes(phy)[node.in, 1]
  	names(originalNodeTimes) <- unique(phy$edge[ ,1])
	internal <- phy$edge[, 1]
	external <- phy$edge[, 2]

	if(length(node) > 0) {

		for (i in 1:length(node)) {
			branchOne <- which(internal == node[i])
			countDes <- ext <- external[branchOne]
			isTerminal <- match(internal, ext)
			terminalCheck <- which(is.na(isTerminal) == FALSE)
			loopStop <- any(terminalCheck)
			
			while (loopStop == TRUE) {
					ext <- external[terminalCheck]
					countDes <- c(countDes, ext)
					isTerminal <- match(internal, ext)
					terminalCheck <- which(is.na(isTerminal) == FALSE)
					loopStop <- any(terminalCheck)
				}

			alreadyCounted <- which(is.na(match(countDes, countDesAll)) == FALSE)
			if (sum(alreadyCounted) > 0) countDes <- countDes[-alreadyCounted]
			countDesAll <- c(countDes, countDesAll)
			}

	tre2 <- phy
	tipsInSample <- countDesAll[which(countDesAll <= Ntip(phy))]
	desInMat <- match(countDesAll, phy$edge[,2])
	tre2$edge <- phy$edge[-desInMat, ]
	keepN <- match(tre2$edge[, 2], node)
	nodeToTips <- tre2$edge[complete.cases(keepN), 2]
	nodalPosition <- which(is.na(keepN) != TRUE)
	originalNodes <- unique(tre2$edge[, 1])
	for (k in 1:length(node)) tre2$edge[which(tre2$edge[ ,2] == node[k]), 2] <- NA
		
	naMat <- which(is.na(tre2$edge[, 2]))
	tipsRemain <- which(tre2$edge[, 2] <= Ntip(phy))
	keeperTips <- tre2$edge[tipsRemain, 2]

	
	tre2$edge[which(tre2$edge[, 2] <= Ntip(phy)), 2] <- NA
	tre2$edge.length <- phy$edge.length[-desInMat]
	internals <- tre2$edge[ ,1]
	tipsInTre2 <- length(which(is.na(tre2$edge[ , 2])))
	n_node <- tipsInTre2 + 1
	nodeTodal <- seq(tre2$edge[1,1], tre2$edge[1,1] + dim(tre2$edge)[1] / 2 - 1, by=1)
	tre2NodeTotal <- match(tre2$edge[,1], nodeTodal)
	nodeTodalNode <- match(nodeTodal, tre2$edge[,1])
	notFoundNode <- nodeTodal[which(is.na(nodeTodalNode))]
	uniqueTre2Int <- unique(tre2$edge[which(is.na(tre2NodeTotal)), 1])
	for(y in 1:length(uniqueTre2Int)) tre2$edge[which(tre2$edge[,1] == uniqueTre2Int[y]), 1] <- notFoundNode[y]
	
	takeValue <- tre2$edge[1, 1] - n_node
	tre2$edge[ ,1] <- tre2$edge[,1] - takeValue
	addOnNodes <- dim(tre2$edge)[1] / 2 - 1
	internalNewTree <- rep(seq(n_node, n_node + addOnNodes, 1), each=2)
	externalNew <- rep(unique(tre2$edge[ ,1]), each=2)
	internalNewMatch <- match(tre2$edge[ ,1], externalNew)
	tre2$edge[,1] <- internalNewTree[internalNewMatch]
	buildInt <- which(complete.cases(tre2$edge[ ,2]) == T)
	extNode <- n_node + 1
	desNodes <- seq(from=extNode, to=extNode + length(buildInt) - 1, 1)
	tre2$edge[buildInt, 2] <- desNodes
	tre2$edge[which(is.na(tre2$edge[,2])), 2] <- 1:tipsInTre2
	apeTree <- list(edge=tre2$edge, tip.label=c(1:tipsInTre2), Nnode=tipsInTre2-1, edge.length=tre2$edge.length)
	class(apeTree) <- "phylo"
	
	toKeep <- apeTree$edge[naMat, 2]
	tipsToKeep <- match(c(1:tipsInTre2), toKeep)
	
	new.tip <- 1:Ntip(apeTree)
	new.tip[-apeTree$edge[naMat, 2]] <- phy$tip.label[-tipsInSample]
	new.tip[apeTree$edge[naMat, 2]] <- paste0("node_", unique(phy$edge[match(nodeToTips, phy$edge[,2]), 2]))
	apeTree$tip.label <- new.tip
	
	if(traits == TRUE) {
		
		keepTipPhy <- match(keeperTips, phy$edge[,2])
		keepTipPhyValue <- nodeEstimate[keepTipPhy]
		newTips <- match(nodeToTips, phy$edge[,2])
		tipObject <- rep(NA, Ntip(apeTree))
		tipObject[which(is.na(tipsToKeep))] <- as.numeric(keepTipPhyValue)	
		needNewTips <- which(is.na(tipObject))
		ageOld <- match(nodeToTips, as.numeric(names(originalNodeTimes)))

		for(p in 1:length(needNewTips)){
			origNode_n <- originalNodeTimes[ageOld][p]		
			nodeN <- nodeToTips[p]
			pastNode <- phy$edge[which(phy$edge[,2] == nodeN), 1]
			pastLength <- phy$edge.length[which(phy$edge[,2] == nodeN)]
			ancAge <- origNode_n + pastLength			
			decAge <- origNode_n
			findNear <- which(c(ancAge - timeCut, - (decAge - timeCut)) == min(c(ancAge - timeCut, - (decAge - timeCut))))
			nearestNode <- which(phy$edge[,2] == c(pastNode, nodeN)[findNear])
			tipObject[needNewTips[p]] <- as.numeric(nodeEstimate[nearestNode])[1]
			}
		
		names(tipObject) <- NULL
		}
	
	} else {
	
	apeTree <- phy	
	tipObject <- nodeEstimate[which(phy$edge[,2] <= Ntip(phy))]
	
	}
	
	node.in.ape <- match(unique(apeTree$edge[,1]), apeTree$edge[,1])
	all.ages.ape <- nodeTimes(apeTree)
	nodeTimesApe <- all.ages.ape[node.in.ape,1]
  	names(nodeTimesApe) <- unique(apeTree$edge[,1])
	tipTimesApe <- all.ages.ape[which(apeTree$edge[,2] <= Ntip(apeTree)), 2] 
	names(tipTimesApe) <- apeTree$edge[which(apeTree$edge[,2] <= Ntip(apeTree)), 2]
	cutOffAge <- nodeTimes(phy)[1,1] - timeCut
	cutOff <- nodeTimes(apeTree)[1,1] - cutOffAge
	preCutOff <- which(tipTimesApe < cutOff)
	preCutOffTips <- match(as.numeric(names(preCutOff)), apeTree$edge[ ,2])
	tipsPreCut <- tipTimesApe[preCutOff] - cutOff
	tipsInApe <- match(preCutOff, apeTree$edge[,2])
	
	apeTree$edge.length[preCutOffTips] <- apeTree$edge.length[preCutOffTips] + tipsPreCut

	if(traits == TRUE) return(list(phy=apeTree, tipData=tipObject))
	if(traits == FALSE) return(phy=apeTree)
}