#' @title Tree plotting for rates
#' @description Plots trees with colours based on rates of trait evolution. Also provides simple coloured plotting for trait values using the \code{\link[ape]{ace}} function in the \pkg{ape} library.
#' @param phy An object of class \code{phylo} (see \pkg{ape})
#' @param x A matrix of trait values.
#' @param traitMedusaObject Output from traitMedusaSummary.
#' @param reconType Colour branches according to rate shifts ("rates" - requires traitMedusaObject) or ancestral state reconstruction ("picReconstruction"  - requires x).
#' @param palette Defines the colour scheme with four options: hotspot.colors (red to blue), heat.colors (yellow to red), cool.colors (blues), combi.colors (yellows to reds and blues)
#' @param type a character string specifying the type of phylogeny to be drawn; it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these.
 #' @param use.edge.length a logical indicating whether to use the edge lengths of the phylogeny to draw the branches (the default) or not (if \code{FALSE}). This option has no effect if the object of class \code{"phylo"} has no `edge.length' element.
#' @param show.tip.label a logical indicating whether to show the tip labels on the phylogeny (defaults to \code{TRUE}, i.e. the labels are shown).
#' @param show.node.label a logical indicating whether to show the node labels on the phylogeny (defaults to \code{FALSE}, i.e. the labels are not shown).
#' @param edge.color a vector of mode character giving the colours used to draw the branches of the plotted phylogeny. These are taken to be in the same order than the component \code{edge} of \code{phy}. If fewer colours are given than the length of \code{edge}, then the colours are recycled.
#' @param edge.width a numeric vector giving the width of the branches of the plotted phylogeny. These are taken to be in the same order than the component \code{edge} of \code{phy}. If fewer widths are given than the length of \code{edge}, then these are recycled.
#' @param edge.lty same than the previous argument but for line types; 1: plain, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash.
#' @param font an integer specifying the type of font for the labels: 1 (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param cex a numeric value giving the factor scaling of the tip and node labels (Character EXpansion). The default is to take the current value from the graphical parameters.
#' @param adj}{a numeric specifying the justification of the text strings of the labels: 0 (left-justification), 0.5 (centering), or 1 (right-justification). This option has no effect if \code{type ="unrooted"}. If \code{NULL} (the default) the value is set with respect of \code{direction} (see details).
#' @param srt a numeric giving how much the labels are rotated in degrees (negative values are allowed resulting in clock-like rotation); the value has an effect respectively to the value of \code{direction} (see Examples). This option has no effect if \code{type = "unrooted"}.
#' @param no.margin a logical. If \code{TRUE}, the margins are set to zero and the plot uses all the space of the device (note that this was the behaviour of \code{plot.phylo} up to version 0.2-1 of `ape' with no way to modify it by the user, at least easily).
#' @param root.edge a logical indicating whether to draw the root edge (defaults to FALSE); this has no effect if `use.edge.length = FALSE' or if `type = "unrooted"'.
#' @param label.offset a numeric giving the space between the nodes and the tips of the phylogeny and their corresponding labels. This option has no effect if \code{type = "unrooted"}.
#' @param underscore a logical specifying whether the underscores in tip labels should be written as spaces (the default) or left as are (if\code{TRUE}).
#' @param x.lim a numeric vector of length one or two giving the limit(s) of the x-axis. If \code{NULL}, this is computed with respect to various parameters such as the string lengths of the labels and the branch lengths. If a single value is given, this is taken as the upper limit.
#' @param y.lim same than above for the y-axis.
#' @param direction a character string specifying the direction of the tree. Four values are possible: "rightwards" (the default),"leftwards", "upwards", and "downwards".
#' @param lab4ut (= labels for unrooted trees) a character string specifying the display of tip labels for unrooted trees: either \code{"horizontal"} where all labels are horizontal (the default), or \code{"axial"} where the labels are displayed in the axis of the corresponding terminal branches. This option has an effect only if \code{type = "unrooted"}.
#' @param tip.color the colours used for the tip labels, eventually recycled (see examples).
#' @return Returns a data frame of colours used in plot along with rate (or ancestral state) range for each colour. 
#' @seealso \code{\link{transformPhylo.ML}}
#' @author Gavin Thomas 
#' @examples
#' # Data and phylogeny
#' data(anolis.tree)
#' data(anolis.data)
#'
#' # female SVL data
#' female.svl <- matrix(anolis.data[,"Female_SVL"], 
#' dimnames=list(rownames(anolis.data)))
#' input.data <- sortTraitData(phy=anolis.tree, y=female.svl, log.trait=TRUE)
#' 
#' # arbitarily reduce data size for speed in this example
#' phy.clade <- extract.clade(input.data[[1]], 182)
#' male.length.clade <- as.matrix(input.data[[2]][match(input.data[[1]]$tip.label, 
#' rownames(input.data[[2]])),])
#' # Identify rate shifts and print and plot results with up to one rate shifts 
#' # and minimum clade size of 10.
#' anolisSVL_MEDUSA <- transformPhylo.ML(male.length.clade, phy=phy.clade, 
#' model="tm1",minCladeSize=10, nSplits=1)
#' anolisSVL_MEDUSA_out <- traitMedusaSummary(anolisSVL_MEDUSA,
#' cutoff=1, AICc=FALSE)
#' colours <- plotPhylo.motmot(phy=phy.clade, traitMedusaObject = anolisSVL_MEDUSA_out,
#' reconType = "rates", type = "fan", cex=0.6, edge.width=3)
#' @import grDevices
#' @export

plotPhylo.motmot <- function(phy, x=NULL, traitMedusaObject=NULL, reconType="rates", type = "phylogram", use.edge.length = TRUE, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0.5, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", palette="hotspot.colors")
    
    
{
   	
	cool.colors <- function (n, alpha = 1) {grDevices::rainbow(n, start = 3/6, end = 4/6, alpha = alpha)}
	heat.colors <- function (n, alpha = 1) {grDevices::rainbow(n, start = 0, end = 1/6, alpha = alpha)}
	hotspot.colors <- function (n, alpha = 1) {grDevices::rainbow(n, start = 0, end = 4/6, alpha = alpha)}  
 	
switch (reconType,
		
	"picReconstruction" = {
		if (is.vector(x)==FALSE) {stop("trait must be a vector with taxon names") }
   	
		a <- ace(x, phy, method="pic")
		reconTrait <- rep(NA, length(phy$edge[,2]))

		# internal nodes   
		reconTrait[match(names(a$ace)[2:length(a$ace)], phy$edge[,2])] <- a$ace[2:length(a$ace)]

		# tips
		idx <- match(names(x), phy$tip.label)
		reconTrait[match(idx, phy$edge[,2])] <- x

		ncolours <- Ntip(phy) + Nnode(phy) -1
		edgeColours <- rep(NA, ncolours)
    
		if (palette=="hotspot.colors"){    use.colors <- rev(hotspot.colors(ncolours))}
		if (palette=="heat.colors"){    use.colors <- rev(heat.colors(ncolours))}
		if (palette=="cool.colors"){    use.colors <- rev(cool.colors(ncolours))}
		if (palette=="combi.colors"){    use.colors <- c(cool.colors(ncolours/2), rev(heat.colors(ncolours/2)))}
		
 
		binsize <- (max(reconTrait) - min(reconTrait)) / length(reconTrait)
		color.bin <- matrix(NA, ncol=length(reconTrait), nrow=2)
		rownames(color.bin) <- c("min", "max")
 
		for (i in 1:length(reconTrait)) { 	
			color.bin[1,i] <- min(reconTrait) + (i-1) * binsize
			color.bin[2,i] <- min(reconTrait) + i* binsize
			}
 
		color.bin[2,length(reconTrait)] <- max(reconTrait)	# necessary due to rounding problems
 	
 	
		for (i in 1: ncolours) {
			color.trait <- (reconTrait[i] >= color.bin["min",]) & (reconTrait[i] <= color.bin["max",])  
			edgeColours[i] <- use.colors[which(color.trait==TRUE)]
			}
		},
		
	"rates" = {
		
		
		if(traitMedusaObject$Rates[1,1] == "0") stop("BM model - no shifts to plot")
		cladeRates <- traitMedusaObject$Rates[,3]
		nodes <- traitMedusaObject$Rates[,1]
		rateType <- traitMedusaObject$Rates[,2]
		
		ncolours <- Ntip(phy) + Nnode(phy) -1
		
		nodeBranch <- data.frame(as.numeric(nodes), rateType)
		
		nodes <- nodeBranch[nodeBranch[,"rateType"]=="clade",1]
		branches <- nodeBranch[nodeBranch[,"rateType"]=="branch",1]
		
   	   	cladeMembers <- cladeIdentity(phy, nodes)
		edgeColours <- rep("black", length(phy$edge[,1]))
		
		lgcladeRates <- log(as.numeric(cladeRates))
		lgcladeRatesNodes <- lgcladeRates[nodeBranch[,"rateType"]=="clade"]
		lgcladeRatesBranches <- lgcladeRates[nodeBranch[,"rateType"]=="branch"]
		
		if (palette=="hotspot.colors"){    use.colors <- rev(hotspot.colors(ncolours))}
		if (palette=="heat.colors"){    use.colors <- rev(heat.colors(ncolours))}
		if (palette=="cool.colors"){    use.colors <- rev(cool.colors(ncolours))}
		if (palette=="combi.colors"){    use.colors <- c(cool.colors(ncolours/2), rev(heat.colors(ncolours/2)))}
		
		binsize <- (max(abs(lgcladeRates)) - -max(abs(lgcladeRates))) / length(edgeColours)
		color.bin <- matrix(NA, ncol=length(edgeColours), nrow=2)
		
		rownames(color.bin) <- c("min", "max")
		
		for (i in 1:length(edgeColours)) { 	
			color.bin[1,i] <- -max(abs(lgcladeRates)) + (i-1) * binsize
			color.bin[2,i] <- -max(abs(lgcladeRates)) + i* binsize
			}
		
		color.bin[2,length(edgeColours)] <- max(abs(lgcladeRates))	# necessary due to rounding problems
		color.bin[1,1] <- -max(abs(lgcladeRates))	# necessary due to rounding problems
		
		for (i in 1:length(nodes)) {
			color.trait <- (lgcladeRatesNodes[i] >= color.bin["min",]) & (lgcladeRatesNodes[i] <= color.bin["max",]) 
			edgeColours[cladeMembers[,i]==1] <- use.colors[which(color.trait==TRUE)]
			}
		
		for (i in 1:length(branches)) {
			color.trait <- (lgcladeRatesBranches[i] >= color.bin["min",]) & (lgcladeRatesBranches[i] <= color.bin["max",]) 
			edgeColours[which(phy$edge[,2]==branches[i])] <- use.colors[which(color.trait==TRUE)]
			}
		
		}
		
		)
		
		
	plot.phylo(x=phy, edge.color = edgeColours, edge.width = edge.width, type=type, use.edge.length=use.edge.length, cex=cex, label.offset=label.offset, show.tip.label=show.tip.label, no.margin=no.margin, direction=direction, adj=adj)
	invisible(as.data.frame(t(rbind(use.colors, exp(color.bin)))))
    }
    
   
   