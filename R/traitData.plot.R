#' @title plot a univariate continuous trait data on a phylogeny
#' @description Plots a phylogeny with lines representing the value of a continuous trait
#' @param y A matrix of trait values with taxon names as rownames.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param col.label colour labels for the traits at the tips and in the histogram
#' @param col.tree colour for the edge labels on the tree
#' @param include.hist Logical. Include a histrogram alongside the plot of the tree?
#' @param cex.plot Numeric. The size of labels for the histogram axis labels
#' @return A plot with the trait values shown at the tips, and a histrogram of the trait values
#' @author Mark Puttick
#' @examples
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' traitData.plot(y=male.length, phy)
#' @export


traitData.plot <- function(y, phy, col.label="red", col.tree="black", cex.plot=0.7, include.hist=F) {

	if(include.hist) {
		par(mfrow=c(1,2), mar=c(3,3,3,3), oma=c(1,1,1,1))
		} else {
		par(mfrow=c(1,1), mar=c(4,4,4,4), oma=c(0,0,0,0))	
		}
			
	max.range <- abs(diff(range(y)))
	trait.data.range <- y + max.range
	if(max(trait.data.range) < 0) 
		trait.data.range <- -trait.data.range
	t.data.scale <- trait.data.range / max(trait.data.range)
	plot(ladderize(phy), show.tip.label=F, edge.col=col.tree)
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    x_cord <- lastPP$xx
    y_cord  <- lastPP$yy
    	tip.x <- x_cord[1:Ntip(phy)]  + (range(x_cord)[2] * 0.005)
	tip.y <-y_cord[1:Ntip(phy)]
	t.data.scale <- (max(x_cord) * t.data.scale) * 0.05
	recycle.color <- rep(col.label, Ntip(phy))
	sapply(1:Ntip(phy), function(x) segments(tip.x[x], tip.y[x], tip.x[x] + t.data.scale[x], tip.y[x], col=col.label[x], xpd=T))
	
	if(include.hist) {
		par(mar=c(0.2,1,0.2,1))
			if(is.null(colnames(y))) {
			name.trait <- "trait"
			} else {
			name.trait <- colnames(y)
			}
		hist(y, col=col.label, xaxs="i", yaxt="n", border="white", main="", xlab="", ylab="", cex.axis=cex.plot, las=2, axes=F)	
		axis(1, tick=F, line=-1.5, cex.axis=cex.plot, las=2)
		mtext("frequency", 4, line=0.2, cex=cex.plot)
		mtext(name.trait, 1, line=0.2, cex=cex.plot)
		}
	}
	
	