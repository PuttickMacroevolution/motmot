#' @title Identify shifts in the rate of trait diversification
#' @description Summarises phenotypic rate variation on phylogenies.
#' @param traitMedusaObject Output of a medusa analysis in \code{\link{transformPhylo.ML}}
#' @param cutoff Cutoff value for differences in AIC scores when comparing models. More complex models with an AIC score more than this number of units lower than simpler models are retained (as per runMedusa in the geiger package).
#' @param AICc If true, AICc is used instead of AIC.
#' @param lowerBound Minimum value for parameter estimates.
#' @param upperBound Maximum value for parameter estimates.
#' @details This functions summarises the output of a "medusa" model in transformPhylo.ML (see below). The best overall model is chosen based on AIC (or AICc if AICc=TRUE). The cut-off point for improvement in AIC score between successively more complex models can be defined using cutoff. The default cutoff is 4 but this is somewhat arbitrary and a "good" cut-off may well vary between data sets so it may well be worth exploring different cutoffs.
#' Summarises fits of "medusa" models ("clade" models generated without any a priori assertion of the location of phenotypic diversification rate shifts). It uses the same approach as the runMedusa function in the geiger package (runMedusa tests for shifts in the rate of lineage diversification). The algorithm first fits a constant-rate Brownian model to the data, it then works iteratively through the phylogeny fitting a two-rate model at each node in turn. Each two-rate model is compared to the constant rate model and the best two-rate model is retained. Keeping the location of this rate shift intact, it then repeats the procedure for a three-rate model and so on. The maximum number of rate shifts can be specified a priori using nSplits. Limits can also be applied to the size (species richness) of clades on which to infer new rate shifts using minCladeSize. This can be useful to enable large trees to be handled but should be used cautiously since specifiying a large minimum clade size may result in biologically interesting nested rate shifts being missed. Equally, very small clade sizes may provide poor estimates of rate that may not be informative.
#' @return ModelFit Summary of the best optimal rate shift model.
#' @return Rates Summary of the rate parameters from the best rate shift model.
#' @return optimalTree A phylo object with branch lengths scaled relative to rate.
#' @seealso \code{\link{transformPhylo.ML}}
#' @references Alfaro ME, Santini F, Brock CD, Alamillo H, Dornburg A, Carnevale G, Rabosky D & Harmon LJ. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. PNAS 106, 13410-13414.
#' O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933
#' Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas
#' @examples 
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' phy.clade <- extract.clade(phy, 182)
#' male.length.clade <- as.matrix(male.length[match(phy.clade$tip.label, rownames(male.length)),])
#' tm1 <- transformPhylo.ML(male.length.clade, phy=phy.clade, model="tm1", minCladeSize=10, nSplits=1)
#' tm1_out <- traitMedusaSummary(tm1, cutoff=1)
#' @export


traitMedusaSummary <- function (traitMedusaObject=NULL, cutoff=4, AICc=TRUE, lowerBound=1e-8, upperBound=200) {
		
	y <- traitMedusaObject[[3]]
	phy <- traitMedusaObject[[4]]
	breaks <- numeric()

    if (AICc==TRUE) { 
    		datCol = 6 
    		} else { 
    		datCol = 5 
    		}
    
	bestModel <- traitMedusaObject[[1]][1,]
	rateType <- traitMedusaObject[[1]][2:nrow(traitMedusaObject[[1]]),2]
	
	for (i in 2:dim(traitMedusaObject[[1]])[1]) {
        if ((traitMedusaObject[[1]][i-1, datCol] - traitMedusaObject[[1]][i, datCol]) < cutoff) 
            break
			bestModel <- traitMedusaObject[[1]][i,]
			where.model <- i
			}

	out <- vector(mode="list", length=3)
	names(out) <- c("ModelFit", "Rates", "optimalTree")
	
	
	
	if (bestModel$node==0) { 
		out$optimalTree <- phy 
		out$ModelFit <- bestModel
		out$Rates <- "Single rate"
		} else {
	
		foo <- function(param) {
			ll <- transformPhylo.ll(y, phyClade, model="clade", nodeIDs=SingleNode, cladeRates=param, rateType=whichRateType)$logLikelihood
			return(as.numeric(ll - as.numeric(bestModel["lnL"]) + 1.92))
		}
		
		
		nodeIDs <- traitMedusaObject[[1]][2:where.model, "node"]
		cladeRates <- traitMedusaObject[[1]][where.model, -c(1:6)]
		if(any(is.na(cladeRates))) cladeRates <- cladeRates[-which(is.na(cladeRates))]
		cladeRates <- as.numeric(cladeRates)
		rateType <- traitMedusaObject[[1]][2:where.model, "shiftPos"]
	
		optimalTree <- transformPhylo(phy, model="clade", nodeIDs=sort(nodeIDs), cladeRates=as.numeric(cladeRates), rateType=rateType)

		out$Rates <- matrix(NA, length(nodeIDs), 5, byrow=TRUE)
		out$Rates <- as.data.frame(out$Rates)
		colnames(out$Rates) <- c("node", "shiftPos", "MLRate", "LowerCI", "UpperCI")
		
		out$ModelFit <- bestModel[ ,3:6]
	
		out$optimalTree <- optimalTree

		for (i in 1:length(nodeIDs)) {
			
			SingleNode <- sort(nodeIDs[i])
			whichRateType <- rateType[i]
			phyClade <- transformPhylo(optimalTree, model="clade", nodeIDs=SingleNode, cladeRates=1/cladeRates[i], rateType=whichRateType)
			LCI <- NULL
			UCI <- NULL
			
			if(foo(lowerBound) < 0 && foo(cladeRates[i]) != foo(lowerBound)) { 
				LCI <- uniroot(foo, interval = c(1e-100, cladeRates[i]))$root 				
				} else {
				LCI <- NA 
				}
			if(foo(upperBound) < 0 && foo(cladeRates[i]) != foo(upperBound)) {
				UCI <- uniroot(foo, interval = c(cladeRates[i], 1000))$root
				} else {
				UCI <- NA
				}
			out$Rates[i,] <- aaa <- c(nodeIDs[i], rateType[i], cladeRates[i], LCI, UCI)
		}
		
		if (any(is.na(out$Rates[,3:4]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}
}
	return(out)
}