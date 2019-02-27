#' @title Identify shifts in the rate of trait diversification through time
#' @description Summarises phenotypic rate variation on phylogenies through 
#' @param timeSliceObject Output of a timeSlice analysis in \code{\link{transformPhylo.ML}}
#' @param cutoff Value for differences in AIC scores when comparing models. More complex models with an AIC score more than this number of units lower than simpler models are retained (as per runMedusa in the geiger package).
#' @param AICc If \code{TRUE}, AICc is used instead of AIC.
#' @param lowerBound Minimum value for parameter estimates.
#' @param upperBound Maximum value for parameter estimates.
#' @param plot.phylo Logical. If \code{TRUE}, the phylogeny is plotted
#' @param cex.tip The character size of tip labels
#' @param tip.offset the distance tip labels should be printed from pendant edges
#' @param phylo.width The width of edges on the phylogeny
#' @param tip.colour Colour of species' tip labels
#' @param colour.ramp The colours signifying different rates from low (first colour) to high (second colour)
#' @details This functions summarises the output of a "timeSlice" model in \code{\link{transformPhylo.ML}} (see below). The best overall model is chosen based on AIC (or AICc if AICc=TRUE). The cut-off point for improvement in AIC score between successively more complex models can be defined using cutoff. The default cutoff is 4 but this is somewhat arbitrary and a "good" cut-off may well vary between data sets so it may well be worth exploring different cutoffs.
#' @return ModelFit Summary of the best optimal rate shift model.
#' @return Rates Summary of the rate parameters from the best rate shift model.
#' @return optimalTree A phylo object with branch lengths scaled relative to rate.
#' @references To Add
#' @author Mark Puttick
#' @seealso \code{\link{transformPhylo.ML}}
#' @import grDevices
#' @examples
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' phy.clade <- extract.clade(phy, 182)
#' male.length.clade <- as.matrix(male.length[match(phy.clade$tip.label, 
#' rownames(male.length)),])
#' timeSlice.10.ml <- transformPhylo.ML(y=male.length.clade, phy=phy.clade, model="timeSlice", 
#' splitTime=c(10))
#' outputSummary <- timeSliceSummary(timeSlice.10.ml, cutoff=0.001, cex.tip=0.5, 
#' phylo.width=2, colour.ramp=c("blue", "red"))
#' @export

timeSliceSummary <- function(timeSliceObject, cutoff=4, AICc=TRUE, lowerBound=1e-8, upperBound=1000, plot.phylo=TRUE, cex.tip=1, tip.offset=1, phylo.width=1, tip.colour="grey50", colour.ramp=c("blue", "red")) {

	if(AICc) {
		diff.AICc <- diff(timeSliceObject$timeSlice[,"AICc"])
	} else {
		diff.AICc <- diff(timeSliceObject$timeSlice[,"AIC"])
	}
	best.model <- which(diff.AICc < -cutoff)
	phy <- timeSliceObject$phy
	y <- timeSliceObject$y
	meserr <- timeSliceObject$meserr
	covPIC <- timeSliceObject$covPIC

	if(length(best.model) == 0) {
		best.mod <- 1
		model.best <- "BM"
		} else {
		best.mod <- max(best.model) + 1
		model.best <- paste("split", best.mod - 1)
		}

		return.top <- timeSliceObject$timeSlice[best.mod, ]
		na.test <- which(is.na(return.top))
		if(length(na.test) > 0) {
		model.param <- return.top[-which(is.na(return.top))]
		} else {
		model.param <- return.top
		}
		
	if(best.mod > 1) {
	
		rate.location <- which(unlist(gregexpr("rate", names(model.param))) != -1)
		rates <- model.param[rate.location]
		time.shift.location <- which(unlist(gregexpr("time", names(model.param))) != -1)
		shift.time <- model.param[time.shift.location]
		lnL <- model.param[1]
		phy.time.slice <- transformPhylo(phy, model="timeSlice", timeRates=rates, splitTime=shift.time)
		LCI <- NULL
		UCI <- NULL
		
		for(i in 1:length(rates)) {
			ratesSingle <-	rates[i]
				
			foo <- function(param) {
					pp <- rates
					pp[i] <- param
					ll <- transformPhylo.ll(y=y, phy=phy, model="timeSlice", timeRates=pp, splitTime=shift.time, meserr = meserr, covPIC = covPIC)$logLikelihood
					return(as.numeric(ll - lnL + 1.92))
				}
							
		if(!any(ratesSingle == lowerBound, ratesSingle == upperBound)) {
			if(foo(lowerBound) < 0) { 
				LCI[i] <- uniroot(foo, interval = c(lowerBound, ratesSingle))$root 				
				} else {
				LCI[i] <- NA 
				}
			if(foo(10000) < 0) {
				UCI[i] <- uniroot(foo, interval = c(ratesSingle, upperBound))$root
				} else {
				UCI[i] <- NA
				}
			}
		}	
			
			
		rates.output <- cbind(rates, cbind(LCI, UCI))
		
	} else {
		shift.time <- NA
		phy.time.slice <- phy
		rates <- 1
		}	
					
	start.time <- nodeTimes(phy)[1,1]
	all.times <- start.time - c(start.time, sort(shift.time, T), 0)
	all.times[1] <- -1
	time.poly <- seq(1, length(all.times))
	
	gradientFun <- grDevices::colorRampPalette(c(colour.ramp[1], colour.ramp[2]))
	
	col.in <- gradientFun(length(rates))[rank(rates, ties.method="first")]

	col.in <- paste0(col.in, "75")

	if(plot.phylo) { 
		par(mfrow=c(2,1), mar=c(2,2,2,2))
		plot(phy, edge.col="#00000000", edge.width=phylo.width, tip.col="white", cex=cex.tip, label.offset=tip.offset)
		for(time.x in time.poly) {
			polygon(c(all.times[time.x], all.times[time.x+1], all.times[time.x+1], all.times[time.x]), c(-1, -1, Ntip(phy)+1, Ntip(phy)+1), col=col.in[time.x], border=F)
		}
		par(new=T)
		plot(phy, edge.col="white", edge.width=phylo.width, tip.col=tip.colour, cex=cex.tip, label.offset=tip.offset)
		##
		all.times.temp <- start.time - c(start.time, sort(shift.time, T), 0)
		all.times2 <- c(0, cumsum(diff(all.times.temp) * rates))
		all.times2[1] <- -1
		if(tail(all.times2, 1) != nodeTimes(phy.time.slice)[1,1]) all.times2[length(all.times2)] <- nodeTimes(phy.time.slice)[1,1]
		plot(phy.time.slice, edge.col="#00000000", edge.width=phylo.width, tip.col="white", cex= cex.tip, label.offset=tip.offset)
		for(time.x in time.poly) {
		polygon(c(all.times2[time.x], all.times2[time.x+1], all.times2[time.x+1], all.times2[time.x]), c(-1, -1, Ntip(phy) + 1, Ntip(phy) + 1), col=col.in[time.x], border=F)
		}
		par(new=T)
		plot(phy.time.slice, edge.col="white", edge.width=phylo.width, tip.col=tip.colour, cex=cex.tip, label.offset=tip.offset)
		if(any(rates == lowerBound)) abline(v=all.times2[which(rates == lowerBound)  + 1], col=col.in[1], lwd=2)
	}	
	output <- list()
	output$ModelFit <- model.best
	output$Rates <- model.param
	if(length(best.model) > 0) output$RatesCI <- rates.output
	output$optimalTree <- phy.time.slice
	return(output)
}