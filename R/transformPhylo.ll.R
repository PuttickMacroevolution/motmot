#' @title Log-likelhood for models of trait evolution. 
#' @description Fits likelihood models for various models of continuous character evolution.
#' @param y A matrix of trait values.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param model The model of trait evolution (see details).
#' @param meserr A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses.
#' @param kappa Value of kappa transform.
#' @param lambda Value of lambda transform.
#' @param delta Value of delta transform.
#' @param alpha Value of alpha (OU) transform.
#' @param psi Value of psi transform.
#' @param lambda.sp Speciation rate estimate for the tree.
#' @param branchLabels Branches on which different psi parameters are estimated in the "multipsi" model.
#' @param nodeIDs Integer - ancestral nodes of clades.
#' @param rateType If model="clade", a vector specifying if rate shift occurs in a clade ("clade") or on the single branch leading to a clade ("branch").
#' @param branchRates Numeric vector specifying relative rates for individual branches
#' @param timeRates The rates (from ancient to recent) for the timeSlice model
#' @param splitTime A split time (measured from the present, or most recent species) at which a shift in the rate occurs for the "timeSlice" model
#' @param acdcRate Value of ACDC transform.
#' @param cladeRates Numeric vector specifying telative rates for clades.
#' @param covPIC Logical. For multivariate analyses, allow for co-variance between traits rates (TRUE) or no covariance in trait rates (FALSE). If FALSE, only the trait variances not co-variances are used.
#' @details This function fits likelihood models (see below) for continuous character evolution where the parameter values are set a priori. The function returns the log-likihood and the Brownian variance (or variance covariance matrix).
#' model="bm"- Brownian motion (constant rates random walk)
#' model="kappa" - fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. 
#' model="lambda" - fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged);
#' model="delta" - fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably.
#'model="free" - fits Mooer's et al's free model where each branch has its own rate of trait evolution. This can be a useful exploratory analysis but it is slow due to the number of parameters, particularly for large trees.
#'model="clade" - fits a model where particular clades are a priori hypothesised to have different rates of trait evolution (see O'Meara et al. 2006; Thomas et al. 2006, 2009). Clades are specified using nodeIDs and are defined as the mrca node. Unique rates for each clade are specified using cladeRates. rateType specifies whether the rate shift occurs in the stem clade or on the single branch leading to the clade.
#' model="OU" - fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal.
#' model="psi" - fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010). 
#' @return brownianVariance Brownian variance (or covariance for multiple traits) given the data and phylogeny
#' @return logLikelihood The log-likelihood of the model and data
#' @seealso \code{\link{transformPhylo.ML}},\code{\link{transformPhylo.MCMC}}, \code{\link{transformPhylo.sim}}
#' @references Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.
#' Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
#' Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. 
#' @references Ingram T. 2011. Speciation along a depth gradient in a marine adaptive radiation. Proc. Roy. Soc. B. 278, 613-618.
#' @references Ingram T,  Harrison AD, Mahler L, Castaneda MdR, Glor RE, Herrel A, Stuart YE, and Losos JB. 2016. Comparative tests of the role of dewlap size in Anolis lizard speciation. Proc. Roy. Soc. B. 283, 20162199. #' Mooers AO, Vamosi S, & Schluter D. 1999. Using phylogenies to test macroevolutionary models of trait evolution: sexual selection and speciation in Cranes (Gruinae). American Naturalist 154, 249-259.
#' O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933
#' Pagel M. 1997. Inferring evolutionary processes from phylogenies. Zoologica Scripta 26, 331-348.
#' Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.
#' Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas, Mark Puttick
#' @examples
#' # Data and phylogeny
#' data(anolis.tree)
#' data(anolis.data)
#' 
#' # anolis.data is not matrix and contains missing data so put together matrix of
#' # relevant traits (here female and male snout vent lengths) and remove species 
#' # with missing data from the matrix and phylogeny
#' anolisSVL <- data.matrix(anolis.data)[,c(5,6)]
#' anolisSVL[,1] <- log(anolisSVL[,1])
#' anolisSVL[,2] <- log(anolisSVL[,2])
#' 
#' tree <- drop.tip(anolis.tree, names(attr(na.omit(anolisSVL), "na.action")))
#' anolisSVL <- na.omit(anolisSVL)
#' 
#' # log likelihood of kappa = 0.1 or 1
#' transformPhylo.ll(anolisSVL, phy=tree, model="kappa", kappa=0.1)
#' transformPhylo.ll(anolisSVL, phy=tree, model="kappa", kappa=1)
#' 
#' # log likelihood of lambda = 0.01 or 1
#' transformPhylo.ll(anolisSVL, phy=tree, model="lambda", lambda=0.01)
#' transformPhylo.ll(anolisSVL, phy=tree, model="lambda", lambda=1)
#' 
#' # log likelihood of delta = 1.5 or 1
#' transformPhylo.ll(anolisSVL, phy=tree, model="delta", delta=1.5)
#' transformPhylo.ll(anolisSVL, phy=tree, model="delta", delta=1)
#' 
#' # log likelihood of alpha = 0.001 or 2
#' transformPhylo.ll(anolisSVL, phy=tree, model="OU", alpha=0.001)
#' transformPhylo.ll(anolisSVL, phy=tree, model="OU", alpha=2)
#' 
#' # log likelihood of psi = 0 (gradual) or 1 (speciational)
#' phy.bd <- birthdeath(tree)
#' mu_over_lambda <- phy.bd[[4]][1]
#' lambda_minus_mu <- phy.bd[[4]][2]
#' lambda.sp <- as.numeric(lambda_minus_mu / (1 - mu_over_lambda))
#' transformPhylo.ll(anolisSVL, phy=tree, model="psi", psi=0, lambda.sp=lambda.sp)
#' transformPhylo.ll(anolisSVL, phy=tree, model="psi", psi=1, lambda.sp=lambda.sp)
#' @export

transformPhylo.ll <- function(y=NULL, phy, model=NULL, meserr=NULL, kappa=NULL, lambda=NULL, delta=NULL, alpha=NULL, psi=NULL, lambda.sp = NULL, nodeIDs=NULL, rateType=NULL, branchRates=NULL, cladeRates=NULL, timeRates=NULL, splitTime=NULL, branchLabels = NULL, acdcRate=NULL,  covPIC = TRUE) {
	
	switch(model,		  
		   
		   "bm" = {
		   transformPhy <- transformPhylo(phy=phy, model="bm", meserr = meserr, y=y)
		   },
		   
		   "kappa" = {
		   transformPhy <- transformPhylo(phy=phy, model="kappa", kappa=kappa, nodeIDs=nodeIDs, meserr = meserr, y=y)
		   },
		   
		   "lambda" = {
		   transformPhy <- transformPhylo(phy=phy, model="lambda", lambda=lambda, meserr = meserr, y=y)
		   },
		   
		   "delta" = {
		   transformPhy <- transformPhylo(phy=phy, model="delta", delta=delta, nodeIDs=nodeIDs, meserr = meserr, y=y)
		   },
		   
		   "free" = {
		   transformPhy <- transformPhylo(phy=phy, model="free", branchRates=branchRates, meserr = meserr, y=y)
		   },
		   
		   "clade" = {
		   transformPhy <- transformPhylo(phy=phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates, rateType=rateType, meserr = meserr, y=y)
		   },
		   
		   "OU" = {
		   transformPhy <- transformPhylo(phy=phy, model="OU", alpha=alpha, nodeIDs=nodeIDs, meserr = meserr, y=y)
		   },
		   
		   "psi" = {
		   transformPhy <- transformPhylo(phy=phy, model="psi", psi=psi, meserr = meserr, y=y, lambda.sp=lambda.sp)
		   },
		   
		   "multipsi" = {
        	transformPhy <- transformPhylo(phy = phy, branchLabels = branchLabels, model = "multipsi", psi = psi, meserr = meserr, y = y, lambda.sp = lambda.sp)
		   },
		   
		   "timeSlice" = {
		   transformPhy <- transformPhylo(phy=phy, model="timeSlice", timeRates=timeRates,  splitTime=splitTime, meserr = meserr, y=y)
		   },
		   
		   	"ACDC" = {
		   transformPhy <- transformPhylo(phy=phy, model="ACDC", acdcRate=acdcRate, nodeIDs=nodeIDs, cladeRates=cladeRates, y=y, meserr = meserr)
		   }
		   
		  )
	
	return(likTraitPhylo(y=y, phy=transformPhy, covPIC = covPIC))
}
