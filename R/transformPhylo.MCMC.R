#' @title Bayesian MCMC for models of trait evolution
#' @description Fits Bayesian models for various models of continuous character evolution using a Metropolis-Hastings Markov Chain Monte Carlo (MCMC) approach
#' @param y A matrix of trait values.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param model The model of trait evolution (see details).
#' @param mcmc.iteration Integer - the number of generations for which to run the MCMC chain
#' @param hiddenSpeciation Logical. If TRUE the psi model will include nodes that are on the 'full.phy' but not the tree pruned of trait data
#' @param useMean Logical. Use the branch-based estimates of extinction of mean (TRUE, default) for the "psi" and "multispi" models only applicable if "hiddenSpeciation" = TRUE
#' @param full.phy The full phylogeny containing the species that do not contain trait data so are not included in 'phy'
#' @param burn.in The proportion of the chain (as given by mcmc.iteration) which to discard as 'burn-in'
#' @param lowerBound Minimum value for parameter estimates
#' @param upperBound Maximum value for parameter estimates
#' @param opt.accept.rate Logical. Perform a pre-run optimisation to achieve an acceptance rate close to 0.44?
#' @param acceptance.sd Numeric. The starting standard deviation for the proposal distribution
#' @param opt.prop The proportion of the mcmc.iteration with which to optimise the acceptance rate.
#' @param fine.tune.bound The distance (+/-) from the optimal acceptance rate of 0.44 at which the fine-tune algorithm will stop. Default = 0.05.
#' @param fine.tune.n The number of iterations with which to optimise the acceptance rate. 
#' @param random.start Use a random starting value for the MCMC run (TRUE), or use the environment set.seed() value
#' @details fine.tune.n The method estimates posterior probabilities using a Metropolis-Hastings MCMC approach. To aide convergence, the model will attempt to reach an acceptable proposal ratio (~0.44) when opt.accept.rate=TRUE. These initial fine-tune repititions only save the standard deviation for the truncated normal distribution that is used for the proposal mechanism. The chain is discarded. Posterior probabilites and MCMC diagnostics come from the seperate output chain that commences after this fine-tune procedure. The MCMC model will estimate the posterior probability for the following models. 
#' \itemize{
#' \item {model="kappa"} {fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. Default bounds from ~0 - 1.}
#' \item {model="lambda"} {fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged). Default bounds from ~0 - 1.}
#' \item {model="delta"} {fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably. Default bounds from ~0 - 5.}
#' \item {model="OU"} {fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal. It is often difficult to distinguish among these possibilities. Default bounds from ~0 - 10.}
#' \item {model="psi"} {fits a acceleration-deacceleration model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010).}
#' \item {model="ACDC"} {fits a model to in which rates can exponentially increased or decrease through time (Blomberg et al. 2003). If the upper bound is < 0, the model is equivalent to the 'Early Burst' model of Harmon et al. 2010. Default rate parameter bounds from ln(1e-10) ~ ln(20) divided by the root age.}
#' }
#' @return median The median estimate of the posterior for the parameter 
#' @return 95.HPD The 95 percent Highest Posterior Density for the parameter
#' @return ESS Effective Sample Size for the posterior
#' @return acceptance.rate The ratio for which new proposals were accepted during the MCMC chain
#' @return mcmc.chain Full MCMC chain containing all iterations (including burn-in)
#' @author Mark Puttick, Gavin Thomas
#' @seealso \code{\link{transformPhylo.ML}}, \code{\link{transformPhylo.ll}}, \code{\link{transformPhylo}}
#' @import coda
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
#' ## please note, this model will be need to run for longer to achieve convergence
#' lambda.mcmc <- transformPhylo.MCMC(y=male.length.clade, phy=phy.clade, 
#' model="lambda", mcmc.iteration=100, burn.in=0.1)
#' @export 

transformPhylo.MCMC <- function(y, phy, model, mcmc.iteration=1000, burn.in=0.1, hiddenSpeciation = FALSE, full.phy=NULL, lowerBound = NULL, upperBound = NULL, opt.accept.rate=TRUE, acceptance.sd=NULL, opt.prop=0.25, fine.tune.bound=0.05, fine.tune.n=30, useMean = FALSE, random.start=TRUE) {

if(random.start) set.seed(as.numeric(Sys.time()))

rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

bounds <- matrix(c(1e-08, 1, 1e-08, 1, 1e-08, 5, 1e-08, 20, 0, 1, 1e-08, 10, 1e-10, 100000), 7, 2, byrow = TRUE)
rownames(bounds) <- c("kappa", "lambda", "delta", "alpha", "psi", "rate", "acdcRate")

if(model == "lambda") {

input.value <- runif(1, 0, 1)

if (is.null(lowerBound)) {
lowerBound <- bounds["lambda", 1]
}
if (is.null(upperBound)) {
upperBound <- bounds["lambda", 2]
}

lik.model <- function(pram) {
lambda.phy <- transformPhylo(phy, model="lambda", y=y, lambda=pram)
return(likTraitPhylo(y, lambda.phy)[[2]])
}

prior.uniform <- function(pram) {
lambda.prior <- dunif(pram, lowerBound, upperBound)
return(sum(lambda.prior))
}

if(is.null(acceptance.sd)) {
stn.dev <- suppressWarnings(diff(range(transformPhylo.ML(y, phy, model="lambda")$Lambda, na.rm=T)))
if(stn.dev == 0) stn.dev <- 1
sd.fine.tune <- stn.dev / 2
} else {
stn.dev <- acceptance.sd
sd.fine.tune <- stn.dev / 2
}	
name.param <- c("Lambda")
}


if(model == "delta") {

input.value <- runif(1, 0, 1)

if (is.null(lowerBound)) {
lowerBound <- bounds["delta", 1]
}
if (is.null(upperBound)) {
upperBound <- bounds["delta", 2]
}

lik.model <- function(pram) {
delta.phy <- transformPhylo(phy, model="delta", y=y, delta=pram)
return(likTraitPhylo(y, delta.phy)[[2]])
}

prior.uniform <- function(pram) {
delta.prior <- dunif(pram, lowerBound, upperBound)
return(sum(delta.prior))
}

if(is.null(acceptance.sd)) {
stn.dev <- suppressWarnings(diff(range(transformPhylo.ML(y, phy, model="delta")$Delta, na.rm=T)))
if(stn.dev == 0) stn.dev <- 1
sd.fine.tune <- stn.dev / 2
} else {
stn.dev <- acceptance.sd
sd.fine.tune <- stn.dev / 2
}	

name.param <- c("Delta")
}

if(model == "kappa") {

input.value <- runif(1, 0, 1)

if (is.null(lowerBound)) {
lowerBound <- bounds["kappa", 1]
}
if (is.null(upperBound)) {
upperBound <- bounds["kappa", 2]
}

lik.model <- function(pram) {
kappa.phy <- transformPhylo(phy, model="kappa", y=y, kappa=pram)
return(likTraitPhylo(y, kappa.phy)[[2]])
}

prior.uniform <- function(pram) {
kappa.prior <- dunif(pram, lowerBound, upperBound)
return(sum(kappa.prior))
}

if(is.null(acceptance.sd)) {
stn.dev <- suppressWarnings(diff(range(transformPhylo.ML(y, phy, model="kappa")$Kappa, na.rm=T)))
if(stn.dev == 0) stn.dev <- 1
sd.fine.tune <- stn.dev / 2
} else {
stn.dev <- acceptance.sd
sd.fine.tune <- stn.dev / 2
}	

name.param <- c("Kappa")
}

if(model == "OU") {

input.value <- runif(1, 0, 1)

if (is.null(lowerBound)) {
lowerBound <- bounds["alpha", 1]
}
if (is.null(upperBound)) {
upperBound <- bounds["alpha", 2]
}

lik.model <- function(pram) {
alpha.phy <- transformPhylo(phy, model="OU", y=y, alpha=pram)
return(likTraitPhylo(y, alpha.phy)[[2]])
}

prior.uniform <- function(pram) {
alpha.prior <- dunif(pram, lowerBound, upperBound)
return(sum(alpha.prior))
}

if(is.null(acceptance.sd)) {
stn.dev <- suppressWarnings(diff(range(transformPhylo.ML(y, phy, model="OU")$Alpha, na.rm=T)))
if(stn.dev == 0) stn.dev <- 1
sd.fine.tune <- stn.dev / 2
} else {
stn.dev <- acceptance.sd
sd.fine.tune <- stn.dev / 2
}	

name.param <- c("alpha")
}

if(model == "ACDC") {

input.value <- runif(1, 0, 1)

rootBranchingTime <- nodeTimes(phy)[1,1]
if (is.null(lowerBound)) {
lowerBound <- log(bounds["acdcRate", 1]) / rootBranchingTime
}
if (is.null(upperBound)) {
upperBound <- log(bounds["acdcRate", 2]) / rootBranchingTime
}

lik.model <- function(pram) {
acdc.phy <- transformPhylo(phy, model="ACDC", y=y, acdcRate=pram)
return(likTraitPhylo(y, acdc.phy)[[2]])
}

prior.uniform <- function(pram) {
acdc.prior <- dunif(pram, lowerBound, upperBound)
return(sum(acdc.prior))
}

if(is.null(acceptance.sd)) {
stn.dev <- suppressWarnings(diff(range(transformPhylo.ML(y, phy, model="ACDC")$acdc, na.rm=T)))
if(stn.dev == 0) stn.dev <- 1
sd.fine.tune <- stn.dev / 2
} else {
stn.dev <- acceptance.sd
sd.fine.tune <- stn.dev / 2
}	

name.param <- c("ACDC.rate")
}

if(model == "psi") {

input.value <- runif(1, 0, 1)

if (is.null(lowerBound)) {
lowerBound <- bounds["psi", 1]
}
if (is.null(upperBound)) {
upperBound <- bounds["psi", 2]
}

if (hiddenSpeciation) {
if (is.null(full.phy)) stop("please provide a full phylogeny")
full.data.match <- match(full.phy$tip.label, rownames(y))
tips.no.data <- full.phy$tip.label[which(is.na(full.data.match))]
phy <- dropTipPartial(full.phy, tips.no.data)
}
 
phy.bd <- birthdeath(phy)
mu_over_lambda <- phy.bd[[4]][1]
lambda_minus_mu <- phy.bd[[4]][2]
lambda.sp <- as.numeric(lambda_minus_mu / (1 - mu_over_lambda))
mu.ext <- as.numeric(lambda_minus_mu / (1 / mu_over_lambda - 1))
		
if (mu.ext > 0) {
phy <- sampleHiddenSp(phy, lambda.sp = lambda.sp, mu.ext = mu.ext, useMean=useMean)
} else {
phy$hidden.speciation <- NULL
}

lik.model <- function(pram) {
psi.phy <- transformPhylo(phy, model="psi", y=y, psi=pram)
return(likTraitPhylo(y, psi.phy)[[2]])
}

prior.uniform <- function(pram) {
psi.prior <- dunif(pram, lowerBound, upperBound)
return(sum(psi.prior))
}

if(is.null(acceptance.sd)) {
stn.dev <- suppressWarnings(diff(range(transformPhylo.ML(y, phy, model="psi")$psi, na.rm=T)))
if(stn.dev == 0) stn.dev <- 1
sd.fine.tune <- 1
} else {
stn.dev <- acceptance.sd
sd.fine.tune <- 1
}	

name.param <- c("psi")
}

model.posterior <- function(pram) return (lik.model(pram) + prior.uniform(pram))	

motmot.mcmc <- function(input.value, iterations, stn.dev, silent=FALSE) {

propose.mcmc <- function(pram) {
return(rtnorm(1, pram, sd=stn.dev, lowerBound, upperBound))
}

mcmc.chain <- matrix(input.value, nrow=1)
    	for (i in 1:iterations) {
        	proposed.move <- propose.mcmc(mcmc.chain[i,])
        	chain.prob <- exp(model.posterior(proposed.move) - model.posterior(mcmc.chain[i,]))
        	if (runif(1) < chain.prob) {
        	mcmc.chain <- rbind(mcmc.chain, proposed.move)
        	} else {
        	mcmc.chain <- rbind(mcmc.chain, mcmc.chain[i,])
        	}
if(!silent) {
cat("\r", "MCMC progress:", sprintf("%.4f", i/iterations * 100), "%")
}
    	}
    	return(mcmc.chain)
}

if(opt.accept.rate) {
cat("optimising acceptance ratio fine-tune")
cat("\n", " ")
cat("running")
count <- 1
stn.dev.2 <- stn.dev
old.ratio <- 10	
best.current <- 0
opt.mcmc <- mcmc.iteration * opt.prop
opt.mcmc.burn <- ceiling(opt.mcmc * burn.in)

opt.done <- FALSE
while(opt.done == FALSE) {
chain <- motmot.mcmc(input.value, opt.mcmc, stn.dev.2, silent=TRUE)
acceptance <- 1 - mean(duplicated(chain[-(1:opt.mcmc.burn),]))	
diff.to.accept <- acceptance - 0.44
new.ratio <- abs(acceptance - 0.44)
cat("\r", "acceptance attempt", signif(acceptance, 3), "best acceptance", signif(best.current, 3), "best SD", signif(stn.dev, 3))
if(new.ratio < old.ratio) {
stn.dev <- stn.dev.2
old.ratio <- new.ratio
best.current <- acceptance
}	
if(abs(acceptance - 0.44) < fine.tune.bound ) {
opt.done <- TRUE
cat("\n", "finished fine.tune")
}
count <- count + 1
if(count > fine.tune.n) {
cat("\n", "finished fine.tune")
opt.done <- TRUE
}
stn.dev.2 <- rtnorm(1, stn.dev, sd.fine.tune, lowerBound, upperBound)
}	
}	
cat("\n", " ")	
chain <- motmot.mcmc(input.value, mcmc.iteration, stn.dev=stn.dev)
burnIn <- ceiling(mcmc.iteration * burn.in)
acceptance.1 <- 1 - mean(duplicated(chain[-(1:burnIn), 1]))
post.burn.in <- chain[-c(1:burnIn), ]

if(dim(chain)[2] == 1) {
ess.mcmc <- coda::effectiveSize(post.burn.in)
median.mcmc <- median(post.burn.in)
hpd.mcmc <- coda::HPDinterval(as.mcmc(post.burn.in))
hpd.mcmc <- as.numeric(hpd.mcmc)
names(ess.mcmc) <- names(median.mcmc) <- name.param
names(hpd.mcmc) <- c("lower 95% HPD", "upper 95% HPD")
} else {
ess.mcmc <- apply(post.burn.in, 2, coda::effectiveSize)
median.mcmc <- apply(chain[-c(1:burnIn), ], 2, median)
hpd.mcmc <- apply(chain[-c(1:burnIn), ], 2, function(x) {
	class(x) <- "mcmc"
	as.numeric(coda::HPDinterval(x))
	})
names(ess.mcmc) <- names(median.mcmc) <- name.param
colnames(hpd.mcmc) <- name.param
}

cat("\n")
output.mcmc <- list(median.mcmc, hpd.mcmc, ess.mcmc, acceptance.1, chain)
rownames(chain) <- NULL
names(output.mcmc) <- c("median", "95.HPD", "ESS", "acceptance.rate", "mcmc.chain")
print(output.mcmc[1:4])
invisible(return(output.mcmc))
}