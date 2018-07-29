#' @title Slice tree (internal function)
#' @description Split tree for the time slice function 
#' @param phy Phylogenetic tree in phylo format
#' @param splitTime Split time in the past
#' @return Branch lengths of sliced tree
#' @author Mark Puttick

sliceTree <- function(phy, splitTime) {
  branchLengths <- c()
  node.times.2 <- node.times <- nodeTimes(phy)
  startTime <- node.times[1, 1]
  split.these.times <- c(splitTime, startTime)
  treeBrLength <- phy$edge.length
  
  for (p in 1:length(split.these.times)) {
    startTime.2 <- startTime - split.these.times[p]
    br.l <- rep(0, dim(node.times)[1])
    these.br <- which(node.times[, 1] > startTime.2)
    br.l[these.br] <- node.times[these.br, 1] - startTime.2
    treeBrLength <- node.times[, 1] - br.l
    intNode <- which(node.times[these.br, 2] != 0)
    test.finish <-
      which(treeBrLength[these.br][intNode] <= node.times[these.br, 2][intNode])
    node.times[these.br, 1] <- treeBrLength[these.br]
    if (length(test.finish) > 0) {
      br.l[these.br][intNode][test.finish] <-
        node.times.2[these.br, 1][intNode][test.finish] - node.times[these.br, 2][intNode][test.finish]
      node.times[these.br, 1][intNode][test.finish] <- 0
    }
    node.times.2 <- node.times
    branchLengths <- cbind(branchLengths, br.l)
  }
  return(branchLengths)
}
