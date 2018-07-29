#' @title Maximum likelihood for models of trait evoluion
#' @description Fits likelihood models for various models of continuous character evolution. Model fitting is based on maximum-likelihood evaluation using phylogenetically independent contrasts. This is exactly equivalent to, but substantially faster than, GLS approaches.
#' @param y A matrix of trait values.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param model The model of trait evolution (see details).
#' @param modelCIs Logical - estimate confidence intervals for parameter estimates.
#' @param nodeIDs Integer - ancestral nodes of clades applicable to rate heterogenous and nested models of evolution (see details)
#' @param rateType If model="clade", a vector specifying if rate shift occurs in a clade ("clade") or on the single branch leading to a clade ("branch").  
#' @param minCladeSize Integer - minimum clade size for inferred rate shift (where model="medusa").
#' @param nSplits Integer - number of rate shifts to apply for model="medusa" and "timeSlice".
#' @param splitTime A split time (measured from the present, or most recent species) at which a shift in the rate occurs for the "timeSlice" model. If splitTime = NULL, then all ages between 1 million year intervals from the root age - 10 Ma to the present + 10 Ma will be included in the search. The best model will be retained in each search, and will be used as a fixed age in the next search. The model will calculate the likelihood for the number of shifts defined by 'nSplits'.
#' @param boundaryAge Only applicable if splitTime = NULL, the age distance from the tips and and youngest tip for which to search for rate shifts. For example, if boundaryAge = 10, only ages between the root age - 10 and the latest tip + 10 will be included in the search. Set to zero to allow testing of all ages.
#' @param testAge If splitTime = NULL, the interval between ages to be tested. For example, if testAge = 1, all 1 Ma ages between the ages defined by 'boundaryAge' will be tested.
#' @param restrictNode List defining monophyletic groups within which no further rate shifts are searched.
#' @param lambdaEst Logical.  Estimate lambda alongside parameter estimates to reduce data noise. Only applicable for models "kappa", "delta", "OU", "psi", "multispi", and "ACDC". Default=FALSE.
#' @param acdcScalar Logical.  For nested EB rate model, simultaneously estimated a rate scalar alongside EB model. Default=FALSE.
#' @param branchLabels Branches on which different psi parameters are estimated in the "multipsi" model
#' @param hiddenSpeciation Logical. If TRUE the psi model will include nodes that are on the 'full.phy' but not the tree pruned of trait data
#' @param full.phy The full phylogeny containing the species that do not contain trait data so are not included in 'phy'
#' @param useMean Logical. Use the branch-based estimates of extinction of mean (TRUE, default) for the "psi" and "multispi" models. Only applicable if "hiddenSpeciation" = TRUE. If FALSE, this will generate a single realisation of the numbers of hidden speciation events on each branch
#' @param profilePlot Logical.  For the single parameter models "kappa", "lambda", "delta", "OU", "psi", "multipsi", and "ACDC", plot the profile of the likelihood.
#' @param lowerBound Minimum value for parameter estimates
#' @param upperBound Maximum value for parameter estimates
#' @param returnPhy Logical. In TRUE the phylogeny with branch lengths transformed by the ML model parameters is returned
#' @param tol Tolerance (minimum branch length) to exclude branches from trait MEDUSA search. Primarily intended to prevent inference of rate shifts at randomly resolved polytomies.
#' @param covPIC Logical. For multivariate analyses, allow for co-variance between traits rates (TRUE) or no covariance in trait rates (FALSE). If FALSE, only the trait variances not co-variances are used.
#' @param meserr A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses. Largely untested - please use cautiously
#' @param n.cores Integer. Set number of computing cores when running model="traitMedusa" (tm1 and tm2 models)
#' @param controlList List. Specify fine-tune parameters for the optim likelihood search
#' @param print.warnings Logical. If TRUE, warnings are issued if confidence intervals fall outside upper or lower bounds
#' @details This function finds the maximum likelihood parameter values for continuous character evolution. For "kappa", "delta", "OU", "multipsi", and "ACDC" it is possible to fit a 'nested' model of evolution in which the ancestral rate of BM swicthes to a different node, as specified by nodeIDs or branchLabels for multipsi. The function returns the maximum-likelihood parameter estimates for the following models.
#' \itemize{
#' \item {model="bm"} {Brownian motion (constant rates random walk).}
#' \item {model="kappa"} {fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. Default bounds from ~0 - 1.}
#' \item {model="lambda"} {fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged). Default bounds from ~0 - 1.}
#' \item {model="delta"} {fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably. If a nodeIDs is supplied, the model will fit a delta model nested within a clade, with a BM fit to the rest of the tree. Default bounds from ~0 - 5.}
#' \item {model="OU"} {fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal. It is often difficult to distinguish among these possibilities. If a nodeIDs is supplied, the model will fit a OU model nested within a clade, with a BM fit to the rest of the tree. For OU models, alternative optimisation are performed with different starting values (1e-8, 0.01, 0.1, 1, 5). Default bounds from ~0 - 10.}
#' \item {model="ACDC"} {fits a model to in which rates can exponentially increased or decrease through time (Blomberg et al. 2003). If the upper bound is < 0, the model is equivalent to the 'Early Burst' model of Harmon et al. 2010. If a nodeIDs is supplied, the model will fit a ACDC model nested within a clade, with a BM fit to the rest of the tree. Default rate parameter bounds from ln(1e-10) ~ ln(20) divided by the root age. Note this process starts on the stem branch leading to the MRCA of the common node, unlike the other methods that start at the common node.}
#' \item {model="psi"} {fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010). Note that the algorithm will automatically estimate speciation and extinction estimates, and will incorporate estimates of 'hidden' speciation if death estimates are greater than 0. }
#' \item {model="multipsi"} {fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate but allows seperate values of psi fitted to seperate branches (Ingram 2010; Ingram et al. 2016). Note that the algorithm will automatically estimate speciation and extinction estimates, and will incorporate estimates of 'hidden' speciation if death estimates are greater than 0.}
#' \item {model="free"} {fits Mooers et al's free model where each branch has its own rate of trait evolution. This can be a useful exploratory analysis but it is slow due to the number of parameters, particularly for large trees. Default rate parameter bounds from ~0 - 200.}
#' \item {model="clade"} {fits a model where particular clades are a priori hypothesised to have different rates of trait evolution (see O'Meara et al. 2006; Thomas et al. 2006, 2009). Clades are specified using nodeIDs and are defined as the mrca node. Default rate parameter bounds from ~0 - 200.}
#' \item {model="tm1"} {fits "clade" models without any a priori assertion of the location of phenotypic diversification rate shifts. It uses the same AIC approach as the runMedusa function in the geiger package (runMedusa tests for shifts in the rate of lineage diversification). The algorithm first fits a constant-rate Brownian model to the data, it then works iteratively through the phylogeny fitting a two-rate model at each node in turn. Each two-rate model is compared to the constant rate model and the best two-rate model is retained. Keeping the location of this rate shift intact, it then repeats the procedure for a three-rate model and so on. The maximum number of rate shifts can be specified a priori using nSplits. Limits can be applied to the size (species richness) of clades on which to infer new rate shifts using minCladeSize. This can be useful to enable large trees to be handled but should be used cautiously since specifiying a large minimum clade size may result in biologically interesting nested rate shifts being missed. Equally, very small clade sizes may provide poor estimates of rate that may not be informative. Limits on the search can also be placed using restrictNode. This requires a list where each element of the list is a vector of tip names that define monophyletic groups. Rate shifts will not be searched for within any of the defined groups. Default rate parameter bounds from ~0 - 1000.}
#' \item {model="tm2"} {this model is similar to "tm1", however, at each node it assesses the fit of two models. The first model is exactly as per "tm1". The second model infers a rate shift on the single branch descending directly from a node but not on any of the descending branches thereafter. Only the best fitting single-branch or whole clade model is retained for the next iteration. If a single-branch shift is favoured, this infers either that there was a rapid shift in trait value along the stem leading to the crown group, or that the members of the clade have undergone parallel shifts. In either case, this can be considered as a change in mean, though separating a single early shift from a clade-parallel shift is not possible with this method. }
#' \item {model="timeSlice"} {A model in which all branch rates change at a time or times set a priori by the user. If  Default rate parameter bounds from ~0 - 1000. If splitTime=NULL, all 1 Ma (as defined by test Age) intervals from the root of the tree - 10 and the youngest tip + 10 will be included in the search. The +/- 10 Ma age can be modified using the argument boundaryAge. At each stage the best fitting model will be stored, and the search will continue until n shifts, with n shifts defined by nSplits. If a single value or vector is used for splitTime, only these ages are included in the search.}
#' }
#' @return Returns the maximum log-likelihood and parameter estimates (with 95 percent confidence intervals if specified). If model="bm" instead returns the Brownian (co)variance and log-likelihood. 
#' @return traitMedusaObject A list in which the first element contains a matrix summarising the parameter estimates and node ids, log-likelihoods, number of parameters (k), AIC and AICc for the best one-rate model, two-rate model, three rate model and so on. The second element is a sub-list where the first element contains all two-rate models, the second element contains all three-rate models and so on. This can be summarised using traitMedusaSummary. The third element is the input trait data. The fourth element is the input phylogeny.
#' @note Confidence intervals are based on the assumption of an asymptotic Chi-square distribution. For multi-parameter models (e.g. rate shift models with more than two rates) the confidence intervals are approximate and are calculated for each parameter in turn while holding all other parameters at their maximum likelihood value. 
#' @seealso \code{\link{transformPhylo.MCMC}}, \code{\link{transformPhylo}}, \code{\link{transformPhylo.ll}}, \code{\link{blomberg.k}}
#' @references Alfaro ME, Santini F, Brock CD, Alamillo H, Dornburg A, Carnevale G, Rabosky D & Harmon LJ. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. PNAS 106, 13410-13414.
#' @references Blomberg SP, Garland T & Ives AR 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57, 717-745.
#' @references Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.
#' @references Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
#' @references Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. 
#' @references Harmon LJ et al. 2010. Early bursts of body size and shape evolution are rare in comparative data. Evolution 57, 717-745.
#' @references Ingram T. 2011. Speciation along a depth gradient in a marine adaptive radiation. Proc. Roy. Soc. B. 278, 613-618.
#' @references Ingram T,  Harrison AD, Mahler L, Castaneda MdR, Glor RE, Herrel A, Stuart YE, and Losos JB. 2016. Comparative tests of the role of dewlap size in Anolis lizard speciation. Proc. Roy. Soc. B. 283, 20162199. 
#' @references Mooers AO, Vamosi S, & Schluter D. 1999. Using phylogenies to test macroevolutionary models of trait evolution: sexual selection and speciation in Cranes (Gruinae). American Naturalist 154, 249-259.
#' @references O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933
#' @references Pagel M. 1997. Inferring evolutionary processes from phylogenies. Zoologica Scripta 26, 331-348.
#' @references Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.
#' @references Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas, Mark Puttick
#' @import MASS
#' @import stats
#' @import coda
#' @import ape
#' @import caper
#' @import stats
#' @import utils
#' @import graphics
#' @import mvtnorm
#' @import parallel
#' @import ks
#' @examples
#' # Data and phylogeny
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' phy.clade <- extract.clade(phy, 182)
#' male.length.clade <- as.matrix(male.length[match(phy.clade$tip.label, rownames(male.length)),])
#'
#' # Brownian motion model
#' transformPhylo.ML(male.length.clade , phy=phy.clade, model="bm")
#' 
#' # Delta
#' transformPhylo.ML(male.length.clade , phy=phy.clade, model="delta", upperBound=2)
#' 
#' # The upper confidence interval for kappa is outside the bounds so try increasing
#' # the upper bound
#' 
#' transformPhylo.ML(male.length.clade , phy=phy.clade, model="delta", upperBound=5)
#' 
#' # Test for different rates in different clades - here with 2 hypothesised
#' # unusual rates compared to the background
#'
#' # This fits the non-censored model of O'Meara et al. (2006)
#' phy.clade$node.label[which(phy.clade$node.label == "3")] <- 2
#' transformPhylo.ML(male.length.clade, phy=phy.clade, model="clade", nodeIDs=c(49, 54))
#' 
#' # Identify rate shifts and print and plot results with upto three rate shifts
#' # and minimum clade size of 20.
#' anolisSVL_MEDUSA <- transformPhylo.ML(male.length.clade, phy=phy.clade, model="tm1",
#' minCladeSize=10, nSplits=2)
#' @export

transformPhylo.ML <- function (y, phy, model = NULL, modelCIs = TRUE, nodeIDs = NULL, rateType = NULL, minCladeSize = 1, nSplits = 2, splitTime = NULL, boundaryAge = 10, testAge = 1, restrictNode = NULL, lambdaEst = FALSE, acdcScalar = FALSE,  branchLabels = NULL, hiddenSpeciation = FALSE, full.phy=NULL, useMean = FALSE, profilePlot = FALSE, lowerBound = NULL, upperBound = NULL, covPIC = TRUE, n.cores = 1, tol = NULL, meserr=NULL, controlList = c(fnscale = -1, maxit = 100, factr = 1e-7, pgtol = 0, type = 2, lmm = 5), returnPhy = FALSE, print.warnings=TRUE) {
   
    bounds <- matrix(c(1e-08, 1, 1e-08, 1, 1e-08, 5, 1e-08, 20, 0, 1, 1e-08, 1000, 1e-10, 20), 7, 2, byrow = TRUE)
    rownames(bounds) <-
      c("kappa", "lambda", "delta", "alpha", "psi", "rate", "acdcRate")
      
      
    lower.function.warning <- function() if(print.warnings) warning("Confidence limits fall outside parameter bounds - consider changing lowerBound")
    upper.function.warning <- function() if(print.warnings) warning("Confidence limits fall outside parameter bounds - consider changing upperBound")  
    
    aic.fun <- function(likelihood, k) return(-2 * likelihood + 2 * k)
    aicc.fun <- function(likelihood, k, n) return(-2 * likelihood + 2 * k + ((2 * k * (k + 1)) / (n - k - 1)))
    x <- NULL
    
    if (acdcScalar && !is.null(nodeIDs))
      upperBound <- -1e6
    
    switch(
      model,
      
      bm = {
        phy <- transformPhylo(phy = phy, model = "bm", meserr = meserr, y = y)
        out <- likTraitPhylo(y, phy, covPIC = covPIC)
        out$logLikelihood <- out[[2]]
        out$brownianVariance <- out[[1]]
        out$root.state <- ancState(phy, y)
        out$AIC <- aic.fun(out$logLikelihood, 2)
        out$AICc <- aicc.fun(out$logLikelihood, 2, Ntip(phy))
        if (returnPhy) {
          out$bmPhy <- phy
        }
      },
      
      kappa = {
        if (is.null(nodeIDs)) {
          nodeIDs <- Ntip(phy) + 1
        } else {
          nodeIDs <- nodeIDs
        }
        
        kappa <- 1
        if (lambdaEst)
          kappa[2] <- 1
        if (is.null(lowerBound)) {
          lowerBound <- bounds["kappa", 1]
          if (lambdaEst) lowerBound[2] <- bounds["lambda", 1]
        }
        if (is.null(upperBound)) {
          upperBound <- bounds["kappa", 2]
          if (lambdaEst) upperBound[2] <- bounds["lambda", 2]
        }
        
        var.funkappa <- function(param) {
          if (length(param) != 2) {
            lambda <- 1
          } else {
            lambda <- param[2]
          }
          kappa <- param[1]
          lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
          return(transformPhylo.ll(y = y, phy = lambdaPhy, kappa = kappa, nodeIDs = nodeIDs, model = "kappa", meserr = meserr, covPIC = covPIC)[[2]])
        }
        
        vo <- optim(kappa, var.funkappa, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
		if (lambdaEst) {
          lambda <- vo$par[2]
        } else {
          lambda <- 1
        }
		lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        
        if (modelCIs == TRUE) {
          kappa.fun <- function(param, chiSq = TRUE) {
            ll <- transformPhylo.ll(y = y, phy = lambdaPhy, kappa = param, model = "kappa", meserr = meserr, covPIC = covPIC)$logLikelihood
            if (chiSq) return(ll - vo$value + 1.92) else return(ll)
          }
          
          if (kappa.fun(lowerBound[1]) < 0) {
            LCI <- uniroot(kappa.fun, interval = c(lowerBound[1], vo$par[1]))$root
          } else {
            LCI <- lowerBound[1]
            lower.function.warning()
          }
          if (kappa.fun(upperBound[1]) < 0) {
            UCI <- uniroot(kappa.fun, interval = c(vo$par[1], upperBound[1]))$root
          } else {
            UCI <- upperBound[1]
            upper.function.warning()
          }
        }
        
        if (profilePlot) {
          par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0))
          kappaCurve <- Vectorize(kappa.fun)
          curve(kappaCurve(x, FALSE), from = lowerBound[1], to = upperBound[1], xlab = expression(kappa), ylab = "log-likelihood", las = 1, main = "profile plot", lwd = 2)
          if (modelCIs) {
            abline(v = c(LCI, vo$par[1], UCI), lty = c(3, 2, 3), lwd = 1, col = "#00000090")
          }
        }
        
        out <- list()
        out$MaximumLikelihood <- vo$value[1]
        if (modelCIs) {
          out$Kappa <- matrix(NA, 1, 3, byrow = TRUE)
          out$Kappa[1,] <- c(vo$par[1], LCI, UCI)
          colnames(out$Kappa) <- c("MLKappa", "LowerCI", "UpperCI")
         } else {
          out$Kappa <- matrix(NA, 1, 1, byrow = TRUE)
          out$Kappa[1,] <- c(vo$par[1])
          colnames(out$Kappa) <- c("MLKappa")
        }
        kappaPhy <- transformPhylo(y = y, phy = lambdaPhy, model = "kappa", kappa = vo$par[1], meserr = meserr)
        out$brownianVariance <- likTraitPhylo(y = y, phy=kappaPhy, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(kappaPhy, y)
        param <- 3
        if (lambdaEst) {
          out$lambda <- vo$par[2]
          param <- 4
        }
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param , Ntip(phy))
        if (returnPhy) {
          out$kappaPhy <- kappaPhy
        }
      },
      
      lambda = {
        lambda <- 1
        if (is.null(lowerBound)) {
          lowerBound <- bounds["lambda", 1]
        }
        if (is.null(upperBound)) {
          upperBound <- bounds["lambda", 2]
        }
        
        var.funlambda <- function(lambda) {
          return(transformPhylo.ll(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr, covPIC = covPIC)[[2]])
        }
		vo <- optim(lambda, var.funlambda, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
        out <- list()
        if (modelCIs == TRUE) {
          lambda.fun <- function(param) {
            ll <- transformPhylo.ll(y, phy, model = "lambda", lambda = param, meserr = meserr, covPIC = covPIC)$logLikelihood
            return(ll - vo$value + 1.92)
          }
          if (lambda.fun(lowerBound) < 0) {
            LCI <- uniroot(lambda.fun, interval = c(lowerBound, vo$par))$root
          } else {
            LCI <- lowerBound
            lower.function.warning()
          }
          if (lambda.fun(upperBound) < 0) {
            UCI <- uniroot(lambda.fun, interval = c(vo$par, upperBound))$root
          } else {
            UCI <- upperBound
            upper.function.warning()
          }     
          out$MaximumLikelihood <- vo$value
          out$Lambda <- matrix(c(vo$par, LCI, UCI), 1, 3, byrow = TRUE)
          colnames(out$Lambda) <- c("MLLambda", "LowerCI", "UpperCI")
        } else {
          out$MaximumLikelihood <- vo$value
          out$Lambda <- matrix(vo$par, 1, 1, byrow = TRUE)
        }
        
        if (profilePlot) {
            par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0))
            lambdaCurve <- Vectorize(lambda.fun)
            curve(lambdaCurve(x), from = lowerBound[1], to = upperBound[1], xlab = expression(paste("Pagel's ", lambda)), ylab = "log-likelihood", las = 1, main = "profile plot", lwd = 2)
            if (modelCIs) {
              abline(v = c(LCI, vo$par[1], UCI), lty = c(3, 2, 3), lwd = 2, col = "#00000090")
            }
          }
          
        lambdaPhy <- transformPhylo(y = y, phy = phy, model = "lambda", lambda = vo$par[1], meserr = meserr)
        out$brownianVariance <- likTraitPhylo(y=y, phy=lambdaPhy, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(lambdaPhy, y)
        out$AIC <- aic.fun(out$MaximumLikelihood, 3)
        out$AICc <- aicc.fun(out$MaximumLikelihood, 3, Ntip(lambdaPhy))
      	if (returnPhy) {
          out$lambdaPhy <- lambdaPhy
        }
 	},
      
      delta = {
        if (is.null(nodeIDs))
          nodeIDs <- Ntip(phy) + 1
        else
          nodeIDs <- nodeIDs
        delta <- 1
        if (lambdaEst)
          delta[2] <- 1
        if (is.null(lowerBound)) {
          lowerBound <- bounds["delta", 1]
          if (lambdaEst)
            lowerBound[2] <- bounds["lambda", 1]
        }
        if (is.null(upperBound)) {
          upperBound <- bounds["delta", 2]
          if (lambdaEst)
            upperBound[2] <- bounds["lambda", 2]
        }
        
        var.fundelta <- function(param) {
          if (length(param) != 2)
            lambda <- 1
          else
            lambda <- param[2]
          delta <- param[1]
          lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
          return(transformPhylo.ll(y = y, phy = lambdaPhy, delta = delta, nodeIDs = nodeIDs, model = "delta", meserr = meserr, covPIC = covPIC)[[2]])
        }
 		vo <- optim(delta, var.fundelta, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
        if (lambdaEst)
          lambda <- vo$par[2]
        else
          lambda <- 1
        lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        
        if (modelCIs == TRUE) {
          delta.fun <- function(param, chiSq = TRUE) {
            ll <- transformPhylo.ll(y = y, phy = lambdaPhy, delta = param, model = "delta", nodeIDs = nodeIDs, meserr = meserr, covPIC = covPIC)$logLikelihood
            if (chiSq) return(ll - vo$value + 1.92) else return(ll)
          }
          if (delta.fun(lowerBound[1]) < 0) {
            LCI <- uniroot(delta.fun, interval = c(lowerBound[1], vo$par[1]))$root
          } else {
            LCI <- lowerBound[1]
            lower.function.warning()
            }
          if (delta.fun(upperBound[1]) < 0) {
            UCI <- uniroot(delta.fun, interval = c(vo$par[1], upperBound[1]))$root
         } else {
            UCI <- upperBound[1]
            upper.function.warning()
            }
        }
        
        if (profilePlot) {
          par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0))
          deltaCurve <- Vectorize(delta.fun)
          curve(deltaCurve(x, FALSE), from = lowerBound[1], to = upperBound[1], xlab = expression(delta), ylab = "log-likelihood", las = 1, main = "profile plot", lwd = 2)
          if (modelCIs)
          	abline(v = c(LCI, vo$par[1], UCI), lty = c(3, 2, 3), lwd = 2, col = "#00000090")
        }
        
        out <- list()
        out$MaximumLikelihood <- vo$value[1]
        if (modelCIs) {
          out$Delta <- matrix(c(vo$par[1], LCI, UCI), 1, 3, byrow = TRUE)
          colnames(out$Delta) <- c("MLDelta", "LowerCI", "UpperCI")
        } else {
          out$Delta <- matrix(vo$par[1], 1, 1, byrow = TRUE)
          colnames(out$Delta) <- c("MLDelta")
        }
        
        deltaPhy <- transformPhylo(y = y, phy = lambdaPhy, model = "delta", delta = vo$par[1], meserr = meserr)
        out$brownianVariance <- likTraitPhylo(phy=deltaPhy, y=y, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(deltaPhy, y)
        param <- 3
        if (lambdaEst) {
          out$lambda <- vo$par[2]
          names(out)[4] <- "Lambda"
          param <- 4
        }
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))
        if (returnPhy) {
          out$deltaPhy <- deltaPhy
        }
      },
      
      OU = {
        if (is.null(nodeIDs))
          nodeIDs <- Ntip(phy) + 1
        else
          nodeIDs <- nodeIDs
        alpha <- c(1e-8, 0.01, 0.1, 1, 5)
        if (lambdaEst)
          alpha[6] <- 1
        if (is.null(lowerBound)) {
          lowerBound <- bounds["alpha", 1]
          if (lambdaEst) lowerBound[2] <- bounds["lambda", 1]
        }
        if (is.null(upperBound)) {
          upperBound <- bounds["alpha", 2]
          if (lambdaEst) lowerBound[2] <- bounds["lambda", 2]
        }
        
        var.funOU <- function(param) {
        	if (length(param) != 2)
            lambda <- 1
          else
            lambda <- param[2]
          alpha <- param[1]
          lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
          return(transformPhylo.ll(y = y, phy = phy, alpha = alpha, nodeIDs = nodeIDs, model = "OU", meserr = meserr, covPIC = covPIC)[[2]])
        }
        
        if (lambdaEst) {
          multiStarts <- lapply(alpha[-7], function(x_multi)
              optim(c(x_multi, 1), var.funOU, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList))
        } else {
          multiStarts <- lapply(alpha, function(x_multi)
              optim(x_multi, var.funOU, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList))
        }
        
        paramOU <- unlist(lapply(multiStarts, function(x_par) x_par[1]))
        likOU <- unlist(lapply(multiStarts, function(x_lik) x_lik[2]))
        if (any(sd(likOU) > 1e-3)) print("Warning - different start values produces different OU likelihoods")
        vo <- multiStarts[[which.max(likOU)]]
		if (lambdaEst)
          lambda <- vo$par[2]
        else
          lambda <- 1
        lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        
        if (modelCIs == TRUE) {
          ou.fun <- function(param, chiSq = TRUE) {
            ll <- transformPhylo.ll(y, lambdaPhy, model = "OU", alpha = param, nodeIDs = nodeIDs, meserr = meserr, covPIC = covPIC)$logLikelihood
            if (chiSq) return(ll - vo$value + 1.92) else return(ll)
          	}
          if (ou.fun(lowerBound[1]) < 0) {
            LCI <- uniroot(ou.fun, interval = c(lowerBound[1], vo$par[1]))$root
          } else {
            LCI <- lowerBound[1]
            lower.function.warning()
          }
          if (ou.fun(upperBound)[1] < 0) {
            UCI <- uniroot(ou.fun, interval = c(vo$par[1], upperBound[1]))$root
          } else {
            UCI <- upperBound[1]
            upper.function.warning()
            }            
        }
        
        if (profilePlot) {
          par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0))
          ouCurve <- Vectorize(ou.fun)
          curve(ouCurve(x, FALSE), from = lowerBound[1], to = upperBound[1], xlab = expression(alpha), ylab = "log-likelihood", las = 1,  main = "profile plot", lwd = 2)
          if (modelCIs) {
            abline(v = c(LCI, vo$par[1], UCI), lty = c(3, 2, 3), lwd = 2, col = "#00000090" )
          }
        }
        
        out <- list()
        out$MaximumLikelihood <- vo$value[1]
        if (modelCIs) {
          out$Alpha <- matrix(c(vo$par[1], LCI, UCI), 1, 3, byrow = TRUE)
          colnames(out$Alpha) <- c("MLAlpha", "LowerCI", "UpperCI")
        } else {
          out$Alpha <- matrix(vo$par[1], 1, 1, byrow = TRUE)
          colnames(out$Alpha) <- c("MLAlpha")
        }
        
        ouPhy <- transformPhylo(y = y, phy = lambdaPhy, model = "OU", alpha = vo$par[1], meserr = meserr)
       	out$brownianVariance <- likTraitPhylo(y=y, phy=ouPhy, covPIC = covPIC)$brownianVariance
        names(out) <- c("MaximumLikelihood", "Alpha", "brownianVariance")
        out$root.state <- ancState(y = y, phy = ouPhy)
        param <- 3
        if (lambdaEst) {
          out$lambda <- vo$par[2]
          param <- 4
        }
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))
        if (returnPhy) out$ouPhy <- ouPhy
      },
      
      ACDC = {
        acdcRate <- -0.01
        if (acdcScalar)
          acdcRate[2] <- 1
        if (lambdaEst)
          acdcRate <- c(acdcRate, 1)
        if (is.null(nodeIDs))
          nodeIDs <- Ntip(phy) + 1
        else
          nodeIDs <- nodeIDs
        rootBranchingTime <- nodeTimes(phy)[1, 1]
        if (is.null(lowerBound)) {
          lowerBound <- log(bounds["acdcRate", 1]) / rootBranchingTime
          if (acdcScalar)
            lowerBound[2] <- 1
          if (lambdaEst)
            lowerBound <- c(lowerBound, bounds["lambda", 1])
        }
        if (is.null(upperBound)) {
          upperBound <- log(bounds["acdcRate", 2]) / rootBranchingTime
          if (acdcScalar)
            upperBound[2] <- 20
          if (lambdaEst)
            upperBound <- c(upperBound, bounds["lambda", 2])
        }
        
        if (acdcScalar) {
          var.funACDC <- function(param) {
            if (lambdaEst) lambda <- tail(param, 1) else lambda <- 1
            acdc.est <- param[1]
            scalarRate <- param[2]
            lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
            return(transformPhylo.ll(lambdaPhy, acdcRate = acdc.est, nodeIDs = nodeIDs, model = "ACDC", y = y, cladeRates = scalarRate, meserr = meserr, covPIC = covPIC)[[2]])
          }
        } else {
          var.funACDC <- function(param) {
            if (lambdaEst) lambda <- param[2] else lambda <- 1
            acdc.est <- param[1]
            lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
            return(transformPhylo.ll(lambdaPhy, acdcRate = acdc.est, nodeIDs = nodeIDs, model = "ACDC", y = y, cladeRates = 1, meserr = meserr, covPIC = covPIC)[[2]])
          }
        }
        vo <- optim(acdcRate, var.funACDC, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
        if (lambdaEst) lambda <- tail(vo$par, 1) else lambda <- 1
        lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        
        if (acdcScalar) cladeRateEst = vo$par[2] else cladeRateEst <- 1
        if (modelCIs) {
          ACDC.fun <- function(param, chiSq = TRUE) {
            ll <- transformPhylo.ll(y, lambdaPhy, acdcRate = param, nodeIDs = nodeIDs, model = "ACDC", cladeRates = cladeRateEst, meserr = meserr, covPIC = covPIC)$logLikelihood
            if (chiSq) return(ll - vo$value + 1.92) else return(ll)
          }
          if (ACDC.fun(lowerBound[1]) < 0) {
            LCI <- uniroot(ACDC.fun, interval = c(lowerBound[1], vo$par[1]))$root
          } else {
            LCI <- lowerBound[1]
            lower.function.warning()
          }
          if (ACDC.fun(upperBound[1]) < 0) {
            UCI <- uniroot(ACDC.fun, interval = c(vo$par[1], upperBound[1]))$root
          } else {
            UCI <- upperBound[1]
            upper.function.warning()
            }
        }
        
        if (profilePlot) {
          par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0))
          acdcCurve <- Vectorize(ACDC.fun)
          curve(acdcCurve(x, FALSE), from = lowerBound[1], to = upperBound[1], xlab = "ACDC rate", ylab = "log-likelihood", las = 1, main = "profile plot", lwd = 2)
          if (modelCIs) {
            abline(v = c(LCI, vo$par[1], UCI), lty = c(3, 2, 3), lwd = 2, col = "#00000090")
          }
        }
        
        out <- list()
        out$MaximumLikelihood <- vo$value[1]
        if (modelCIs) {
          out$ACDC <- matrix(c(vo$par[1], LCI, UCI), 1, 3, byrow = TRUE)
          colnames(out$ACDC) <- c("MLacdc", "LowerCI", "UpperCI")
        } else {
          out$ACDC <- matrix(vo$par[1], 1, 1, byrow = TRUE)
          colnames(out$ACDC) <- c("MLacdc")
        }
        if (acdcScalar) out$scalar <- scaleClade <- vo$par[2] else scaleClade <- 1
          
        acdcPhy <- transformPhylo( y = y, phy = lambdaPhy, model = "ACDC", acdcRate = vo$par[1], nodeIDs = nodeIDs, cladeRates = scaleClade, meserr = meserr)  
        out$brownianVariance <- likTraitPhylo(phy=acdcPhy, y=y, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(phy=acdcPhy, y=y)
        param <- length(vo$par) + 2
        if (lambdaEst) {
          out$lambda <- tail(vo$par, 1)
          param <- param + 1
        }
        n <- Ntip(phy)
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))
        if (returnPhy) out$acdcPhy <- acdcPhy
      },
      
      psi = {
       
        psi <- 1
        if (lambdaEst) psi[2] <- 1
        if (is.null(lowerBound)) {
          lowerBound <- bounds["psi", 1]
          if (lambdaEst)
            lowerBound[2] <- bounds["lambda", 1]
        }
        if (is.null(upperBound)) {
          upperBound <- bounds["psi", 2]
          if (lambdaEst)
            upperBound[2] <- bounds["lambda", 2]
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
        
        var.funpsi <- function(param) {
        	if (length(param) != 2) lambda <- 1 else lambda <- param[2]
        psi <- param[1]
        lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        return(transformPhylo.ll(y = y, phy = lambdaPhy, psi = psi, model = "psi", meserr = meserr, covPIC = covPIC, lambda.sp = lambda.sp)[[2]])
        }
        vo <- optim(psi, var.funpsi, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
        
        if (lambdaEst) lambda <- vo$par[2] else lambda <- 1
        lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        
        if (modelCIs == TRUE) {
          psi.fun <- function(param, chiSq = TRUE) {
            ll <- transformPhylo.ll(y, lambdaPhy, model = "psi", psi = param, meserr = meserr, covPIC = covPIC, lambda.sp = lambda.sp)$logLikelihood
            if (chiSq) return(ll - vo$value + 1.92) else return(ll)
			}
          if (psi.fun(lowerBound[1]) < 0) {
            LCI <- uniroot(psi.fun, interval = c(lowerBound[1], vo$par[1]))$root
          } else {
            LCI <- lowerBound[1]
            lower.function.warning()
            }
          if (psi.fun(upperBound[1]) < 0) {
            UCI <- uniroot(psi.fun, interval = c(vo$par[1], upperBound[1]))$root
          } else {
            UCI <- upperBound[1]
            upper.function.warning()
            }
        }
        
        if (profilePlot) {
          par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0))
          psiCurve <- Vectorize(psi.fun)
          curve(psiCurve(x, FALSE), from = lowerBound[1], to = upperBound[1], xlab = expression(psi), ylab = "log-likelihood", las = 1, main = "profile plot", lwd = 2)
          if (modelCIs) {
            abline(v = c(LCI, vo$par[1], UCI), lty = c(3, 2, 3), lwd = 2, col = "#00000090")
          }
        }
        
        out <- list()
        out$MaximumLikelihood <- vo$value[1]
        if (modelCIs) {
          out$psi <- matrix(c(vo$par[1], LCI, UCI), 1, 3, byrow = TRUE)
          colnames(out$psi) <- c("MLpsi", "LowerCI", "UpperCI")
        } else {
          out$psi <- matrix(vo$par[1], 1, 1, byrow = TRUE)
          colnames(out$psi) <- c("MLpsi")
        }
        psiPhy <- transformPhylo(y=y, phy=lambdaPhy, model="psi", psi = vo$par[1], meserr = meserr, lambda.sp = lambda.sp)
        out$brownianVariance <- likTraitPhylo(y=y, phy=psiPhy, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(phy=psiPhy, y=y)
        param <- 3
        if (lambdaEst) {
          out$lambda <- vo$par[2]
          param <- 4
        }
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))
        if (returnPhy) out$psiPhy <- psiPhy
      },
      
     multipsi = {
		if (is.null(branchLabels)) stop("for 'multipsi' model must provide branchLabels giving state for each branch")
		states <- levels(factor(branchLabels))
		start <- rep(0.1, length(states))
		if (lambdaEst) start <- c(start, 1)
        if (is.null(lowerBound)) {
            lowerBound <- bounds[rep("psi",length(states)), 1]
            if (lambdaEst) lowerBound[length(states) + 1] <- bounds["lambda", 1]
        }
        if (is.null(upperBound)) {
            upperBound <- bounds[rep("psi",length(states)), 2]
            if (lambdaEst) upperBound[length(states) + 1] <- bounds["lambda", 2]
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
        
        var.funmultipsi <- function(param) {
        		all.param <- length(param)
			if (lambdaEst) {
				lambda <- param[all.param]
				psi <- param[-all.param]
				} else {
        			psi <- param
        			lambda <- 1
        			}
        			lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
            	return(transformPhylo.ll(y = y, phy = lambdaPhy, branchLabels = branchLabels, psi = psi, model = "multipsi", meserr = meserr, covPIC = covPIC, lambda.sp = lambda.sp)[[2]])
        }
        vo <- optim(start, var.funmultipsi, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
        
        if (lambdaEst) lambda <- tail(vo$par, 1) else lambda <- 1
        lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)

        out <- list()
        out$MaximumLikelihood <- vo$value[1]	
        out$psi <- matrix(NA, length(states), 3, byrow = TRUE,dimnames=list(states,c("MLpsi", "LowerCI", "UpperCI")))
        out$psi[,1] <- vo$par[1:length(states)]        


        if (modelCIs == TRUE) {
            psi.fun <- function(param, chiSq = TRUE) {
				psi <- as.numeric(vo$par)
				psi[i] <- param
                ll <- transformPhylo.ll(y, lambdaPhy, model = "multipsi", branchLabels = branchLabels, psi = psi, meserr = meserr, covPIC = covPIC, lambda.sp = lambda.sp)$logLikelihood
                if(chiSq) return(ll - vo$value + 1.92) else return(ll)
            }        	
        	for (i in 1:length(states)) { 
        		 if (psi.fun(lowerBound[i]) < 0) {
                    LCI <- uniroot(psi.fun, interval = c(lowerBound[i], vo$par[i]))$root
                } else {
                    LCI <- NA
                }
                if (psi.fun(upperBound[i]) < 0) {
                    UCI <- uniroot(psi.fun, interval = c(vo$par[i], upperBound[i]))$root
                } else {
                    UCI <- NA
                }
                out$psi[i, 2:3] <- c(LCI, UCI)
            }
		}
		
		if (profilePlot) {
          par(mar = c(5,5,5,5), oma = c(0, 0, 0, 0), mfrow=c(length(states), length(states)))
          	for (i in 1:length(states)) {
          		psiCurve <- Vectorize(psi.fun)
          		curve(psiCurve(x, FALSE), from = lowerBound[1], to = upperBound[1], xlab = expression(psi), ylab = "log-likelihood", las = 1, main = paste0("profile plot psi ", i), lwd = 2)
          		if (modelCIs) {
            		abline(v = c(out$psi[i,2], out$psi[i,1], out$psi[i,3]), lty = c(3, 2, 3), lwd = 2, col = "#00000090")
          			}
        		} 	
         	}      
		multipsiPhy <- transformPhylo(y=y, phy=lambdaPhy, model="multipsi", psi = vo$par[1:length(states)], meserr = meserr, lambda.sp = lambda.sp, branchLabels = branchLabels)
        out$brownianVariance <- likTraitPhylo(y=y, phy=multipsiPhy, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(phy=multipsiPhy, y=y)
        param <- length(states) + 2
        	if (lambdaEst) {
          		out$lambda <- vo$par[length(states) + 1]
          		param <- param + 1
        		}   
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))
        if (returnPhy) out$psiPhy <- multipsiPhy
	},
      
     free = {
        if (is.null(lowerBound)) {
            lowerBound <- bounds["rate", 1]
        }
        if (is.null(upperBound)) {
            upperBound <- bounds["rate", 2]
        }
        branchRates <- rep(1, length(phy$edge.length))
        var.funfree <- function(branchRates) {
            return(transformPhylo.ll(y, phy, branchRates = branchRates, model = "free", meserr = meserr, covPIC = covPIC)[[2]])
        }
        vo <- optim(branchRates, var.funfree, method = "L-BFGS-B", lower = 0, control = c(fnscale = -1, maxit = 10, factr = 1e+14))
        phy2 <- phy
        phy2$edge.length <- phy$edge.length * vo$par
        out <- vector(mode = "list", length = 3)
        output.par <- likTraitPhylo(y=y, phy=phy2, covPIC = covPIC)
        out$MaximumLikelihood <- output.par$logLikelihood
        out$brownianVariance <- output.par$brownianVariance
        out$Rates <- vo$par
        if (vo$convergence == 0) {
            out$Convergence <- "Successful"
        } else {
            out$Convergence <- "Failed"
        }
        param <- length(branchRates) + 2
        if (lambdaEst) {
          out$lambda <- vo$par[length(branchRates) + 1]
          param <- param + 1
        }
		out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))       
        if (returnPhy) out$freePhy <- phy2
    },
      
      clade = {
        if (is.null(lowerBound)) {
          lowerBound <- bounds["rate", 1]
        }
        if (is.null(upperBound)) {
          upperBound <- bounds["rate", 2]
        }
        
        if (is.null(rateType)) {
          rateType <- rep("clade", length(nodeIDs))
        } else {
          rateType <- rateType
        }
        
        cladeRates <- rep(1, length(nodeIDs))
        var.funclade <- function(param) {
          return(transformPhylo.ll(y, phy, nodeIDs = nodeIDs, cladeRates = param, model = "clade", rateType = rateType)[[2]])
        }
        vo <- optim(cladeRates, var.funclade, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
        phyCladeFull <- transformPhylo(phy, model = "clade", nodeIDs = nodeIDs, rateType = rateType, cladeRates = vo$par, y = y, meserr = meserr)
        out <- list()
        if (modelCIs == TRUE) {
          free.fun <- function(param) {
            ll <- transformPhylo.ll(y, phyClade, model = "clade", nodeIDs = SingleNode, rateType = whichRateType, cladeRates = param,  meserr = meserr, covPIC = covPIC)$logLikelihood
            return(ll - vo$value + 1.92)
          }
          out$Rates <- matrix(vo$value, length(nodeIDs), 4, byrow = TRUE)
          colnames(out$Rates) <- c("node", "MLRate", "LowerCI", "UpperCI")
          for (i in 1:length(nodeIDs)) {
            SingleNode <- nodeIDs[i]
            whichRateType <- rateType[i]
            phyClade <- transformPhylo(phyCladeFull, model = "clade", nodeIDs = SingleNode, rateType = whichRateType, cladeRates = 1 / vo$par[i], y = y, meserr = meserr)
            if (free.fun(lowerBound) < 0) {
              LCI <- uniroot(free.fun, interval = c(lowerBound, vo$par[i]))$root
            } else {
              LCI <- lowerBound
              lower.function.warning() 
            }
            if (free.fun(upperBound) < 0) {
              UCI <- uniroot(free.fun, interval = c(vo$par[i], upperBound))$root
            } else {
              UCI <- upperBound
              upper.function.warning()
            }
            out$Rates[i,] <- c(nodeIDs[i], vo$par[i], LCI, UCI)
          }
        } else {
          out$Rates <- matrix(NA, length(nodeIDs), 2, byrow = TRUE)
          colnames(out$Rates) <- c("node", "MLRate")
          out$Rates[, 1] <- nodeIDs
          out$Rates[, 2] <- vo$par
         }
        out$MaximumLikelihood <- vo$value
        out$brownianVariance <- likTraitPhylo(y = y, phy = phyCladeFull, covPIC = covPIC)$brownianVariance
        out$root.state <- ancState(phyCladeFull, y)
        param <- 2 + length(vo$par)
        out$AIC <- aic.fun(out$MaximumLikelihood, param)
        out$AICc <- aicc.fun(out$MaximumLikelihood, param, Ntip(phy))
        if (returnPhy) out$freePhy <- phyCladeFull  
      },
      
    tm1 = {
			if(is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
			}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
			}
			nodes <- unique(sort(phy$edge[, 2]))
			nodeDepth <- node.depth(phy)
			nodesByRichness <- cbind(richness = nodeDepth[nodes], node = nodes)
			searchNode <- nodesByRichness[nodesByRichness[, 1] >= minCladeSize, 2]
			if (is.null(restrictNode) == FALSE) {
				ngroups <- length(restrictNode)
				cm <- clade.matrix(phy)
				colnames(cm$clade.matrix) <- cm$tip.label
				cmMat <- cm$clade.matrix
				skipNodes <- numeric()
				for (i in 1:ngroups) {
					matDrop <- cmMat[, restrictNode[[i]]]
					nodeDrop <- rowSums(matDrop) < length(restrictNode[[i]]) & rowSums(matDrop) >= 1
					skipNodes <- c(skipNodes, as.numeric(rownames(matDrop)[nodeDrop]))
					}
				skipNodes <- unique(c(skipNodes, phy$edge[phy$edge.length < tol, 2]))
				searchNode <- setdiff(searchNode, skipNodes)
			} else {
				skipNodes <- unique(phy$edge[phy$edge.length < tol, 2])
				searchNode <- setdiff(searchNode, skipNodes)
			}

			n <- length(y)
			BERateOut <- matrix(NA, nrow = nSplits + 1, ncol = (6 + nSplits))
			fullModelOut <- vector(mode = "list", length = nSplits)
			MLsingle <- transformPhylo.ll(y=y, phy=phy, model="bm", meserr = meserr, covPIC = covPIC)[[2]]
			AICsingle <- aic.fun(MLsingle, 2)
			AICcsingle <- aicc.fun(MLsingle, 2, n)
			colnames(BERateOut) <- 1:dim(BERateOut)[2]
			colnames(BERateOut)[1:6] <- c("node", "shiftPos", "lnL", "n.params", "AIC", "AICc")
			colnames(BERateOut)[-c(1:6)] <- paste0("rate.", 1:nSplits)
			rownames(BERateOut) <- c("BM     ", paste0("shift.", 1:nSplits))
			BERateOut[1, 1:6] <- c(0, 1, MLsingle, 2, AICsingle, AICcsingle)			
			BERateOut <- data.frame( BERateOut)
			cat("\n", "BM model")
			cat("\n")
			print(BERateOut[1,1:6])
	
			for (k in 1:nSplits) {
				cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 1]))
				if (k == 1) {
					searchNode <- searchNode
					} else { 
					searchNode <- searchNode[!(searchNode %in% bestNodes)]
					}
				for (i in searchNode) {
					if (k == 1) {
						currentNodeIDs <- i 
						} else {
						currentNodeIDs <- c(bestNodes, i)
						}
					cladeRates <- rep(1, length(currentNodeIDs))
					var.funclade.tm1 <- function(cladeRates) {
					return(transformPhylo.ll(y, phy, nodeIDs = currentNodeIDs, rateType = rep("clade", length(currentNodeIDs)), cladeRates = cladeRates, model = "clade", meserr = meserr, covPIC = covPIC)[[2]])
					}
					currentCladeModel <- optim(cladeRates, var.funclade.tm1, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
					fullModelOut[[k]] <- rbind(fullModelOut[[k]], c(node = as.integer(i), shiftPos = 1, ML = currentCladeModel$value, currentCladeModel$par))
					currentModel <- currentCladeModel
					shiftPos = 1
					param <- (k * 2) + 1
					currentModel <- list(currentModel, i, cladeMembers)
					currentML <- currentModel[[1]]$value
					AIC <- aic.fun(currentModel[[1]]$value,  param)
					AICc <- aicc.fun(currentModel[[1]]$value, param, n)

					if (i == min(searchNode)) {
						BERateOut[k+1, 1:(6 + k)] <- c(currentModel[[2]], shiftPos, currentModel[[1]]$value, param, AIC, AICc, rev(currentModel[[1]]$par))
					} else {
						if (currentML > BERateOut[k+1, 3]){
						BERateOut[k+1, 1:(6 + k)] <- c(currentModel[[2]], shiftPos, currentModel[[1]]$value, param, AIC, AICc, rev(currentModel[[1]]$par))
						}
					}
				}

			if (k == 1) bestNodes <- BERateOut[k + 1, 1] else bestNodes <- c(bestNodes, BERateOut[k + 1, 1])
			if (sum(BERateOut[, "shiftPos"] == 1, na.rm=T) > 0) BERateOut[which(BERateOut[, "shiftPos"] == 1), "shiftPos"] <- "clade"
			cat("\n", "Shift", k)
			cat("\n")
			print(BERateOut[k+1,c(1:(6+k)) ])
			}
			out <- list(as.data.frame(BERateOut), fullModelOut, y, phy)
		}, 
      
      tm2 = {
			if (is.null(lowerBound)) lowerBound <- bounds["rate", 1]
			if (is.null(upperBound))upperBound <- bounds["rate", 2]
			nodes <- unique(sort(phy$edge[, 2]))
			nodeDepth <- node.depth(phy)
			nodesByRichness <- cbind(richness = nodeDepth[nodes], node = nodes)
			searchNode <- nodesByRichness[nodesByRichness[, 1] >=  minCladeSize, 2]
			searchNodeType <- c(searchNode, searchNode[searchNode > Ntip(phy)])
			searchType <- c(rep("clade", length(searchNode)), rep("branch", sum(searchNode > Ntip(phy))))
			searchNodeTypeAll <- data.frame(searchType, searchNodeType, stringsAsFactors = FALSE)
			if (is.null(restrictNode) == FALSE) {
				ngroups <- length(restrictNode)
				cm <- clade.matrix(phy)
				colnames(cm$clade.matrix) <- cm$tip.label
				cmMat <- cm$clade.matrix
				skipNodes <- numeric()
				for (i in 1:ngroups) {
					matDrop <- cmMat[, restrictNode[[i]]]
					nodeDrop <- rowSums(matDrop) < length(restrictNode[[i]]) & rowSums(matDrop) >= 1
					skipNodes <- c(skipNodes, as.numeric(rownames(matDrop)[nodeDrop]))
					}
				skipNodes <- unique(c(skipNodes, phy$edge[phy$edge.length < tol, 2]))
				searchNode <- setdiff(searchNode, skipNodes)
			} else {
				skipNodes <- unique(phy$edge[phy$edge.length < tol, 2])
				searchNode <- setdiff(searchNode, skipNodes)
			}
			
			tm2.fun <- function(whichSearchNode) which(searchNodeTypeAll[, 2] == whichSearchNode)
			searchNodeTypeAll <- searchNodeTypeAll[sort(unlist(sapply(searchNode, tm2.fun))), ]
			n <- length(y)
			BERateOut <- matrix(NA, nrow = nSplits + 1, ncol = (6 + (nSplits)))
			fullModelOut <- vector(mode = "list", length = nSplits)
			bestModel <- matrix(NA, ncol = 2, nrow = nSplits)
			fullModelOut <- vector(mode = "list", length = nSplits)
			MLsingle <- transformPhylo.ll(y=y, phy=phy, model="bm", meserr = meserr, covPIC = covPIC)[[2]]
			AICsingle <- aic.fun(MLsingle, 2)
			AICcsingle <- aicc.fun(MLsingle, 2, Ntip(phy))
			BERateOut[1, 1:6] <- c(0, 1, MLsingle, 2, AICsingle,AICcsingle)
			colnames(BERateOut) <- 1:dim(BERateOut)[2]
			colnames(BERateOut)[1:6] <- c("node", "shiftPos", "lnL", "n.params", "AIC", "AICc")
			colnames(BERateOut)[-c(1:6)] <- paste0("rate.", 1:nSplits)
			rownames(BERateOut) <- c("BM     ", paste0("shift.", 1:nSplits))
			BERateOut[1, 1:6] <- c(0, 1, MLsingle, 2, AICsingle, AICcsingle)			
			BERateOut <- data.frame( BERateOut)
			cat("\n", "BM model")
			cat("\n")
			print(BERateOut[1,1:6])

			for (k in 1:nSplits) {
				cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 1]))
				bestModel[bestModel[, 2] == "1", 2] <- "clade"
				bestModel[bestModel[, 2] == "2", 2] <- "branch"
				if (k == 1) {
					searchNodeTypeAll <- searchNodeTypeAll
				} else {
					dropNodeRow <- which(searchNodeTypeAll[, 1] == bestModel[k - 1, 2] & searchNodeTypeAll[, 2] == as.numeric(bestModel[k - 1, 1]))
					rowDropLogical <- row(searchNodeTypeAll)[, 1] != dropNodeRow
					searchNodeTypeAll <- searchNodeTypeAll[rowDropLogical, ]
				}

### Weird stuff going on here. Updating of currentNodeIDs seems to be the problem.Needs to be taken out of the function? Compare old with new version to see where problem is.

				if (k == 1) currentNodeIDs <- NULL else currentNodeIDs <- na.omit(as.numeric(bestModel[, 1]))
				if (k == 1) currentRateType <- NULL else currentRateType <- na.omit(bestModel[, 2])

				foo2 <- function(x_foo, k, currentNodeIDs=currentNodeIDs, currentRateType=currentRateType) {
					currentNodeIDs <- c(currentNodeIDs, x_foo[,2])
					currentRateType <- c(currentRateType, x_foo[,1])
					cladeRates <- rep(1, length(currentNodeIDs))
					var.funclade.tm2 <- function(cladeRates) return(transformPhylo.ll(y, phy, nodeIDs = currentNodeIDs, rateType = currentRateType, cladeRates = cladeRates, model = "clade", meserr = meserr, covPIC = covPIC)[[2]])
					currentCladeModel <- optim(cladeRates, var.funclade.tm2, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
					if (x_foo[,1] == "clade") shiftPos <- 1
					if (x_foo[,1] == "branch") shiftPos <- 2
					currentModel <- currentCladeModel
					param <- k * 2 + 1
					currentModel <- list(currentModel, x_foo[2], cladeMembers)
					currentML <- currentModel[[1]]$value
					AIC <- aic.fun(currentModel[[1]]$value, param)
					AICc <- aicc.fun(currentModel[[1]]$value, param, n)
					out <- c(currentModel[[2]][,1], shiftPos, currentModel[[1]]$value, param, AIC, AICc, currentModel[[1]]$par)
					return(out)
					}

				searchNodeList <- vector(mode="list", length=length(searchNodeTypeAll[,1]))
				for (i in 1:length(searchNodeTypeAll[,1])) searchNodeList[[i]] <- searchNodeTypeAll[i,]
				all.output <- parallel::mclapply(searchNodeList, FUN=foo2, k=k, currentNodeIDs=currentNodeIDs, currentRateType=currentRateType,  mc.cores=n.cores)
				all.output <- matrix(unlist(all.output), nrow=length(searchNodeList) , byrow=TRUE)
				best.out <- all.output[which.max(all.output[,3]),]
				fullModelOut[[k]] <- all.output
				BERateOut[k+1, 1:(6 + k)] <- best.out
				bestModel[k,] <- unlist(BERateOut[k+1, 1:2])				
				if (sum(BERateOut[, "shiftPos"] == 1, na.rm=T) > 0) BERateOut[which(BERateOut[, "shiftPos"] == 1), "shiftPos"] <- "clade"
				if (sum(BERateOut[, "shiftPos"] == 2, na.rm=T) > 0) BERateOut[which(BERateOut[, "shiftPos"] == 2), "shiftPos"] <- "branch"
				cat("\n", "Shift", k)
				cat("\n")
				print(BERateOut[k+1,c(1:(6+k)) ])
				}
				out <- list(as.data.frame(BERateOut), fullModelOut, y, phy)
		},
      
      timeSlice = {
        maxTime <- nodeTimes(phy)[1, 1]
        n <- Ntip(phy)
        if (is.null(controlList["pgtol"])) controlList["pgtol"] <- 0.1
        if (is.null(lowerBound))
          lowerBound <- bounds["rate", 1]
        if (is.null(upperBound))
          upperBound <- bounds["rate", 2]
        
        if (is.null(splitTime)) {
          splitTime <- seq(maxTime - boundaryAge, boundaryAge,-testAge)
          output.mat <- matrix(NA, nrow = nSplits + 1, ncol = 6 +  (2 * nSplits))
          bm.model <- likTraitPhylo(y, phy, covPIC = covPIC)
          log.lik <- as.numeric(bm.model[[2]])
          aic <- aic.fun(log.lik, 2)
          aicc <- aicc.fun(log.lik, 2, Ntip(phy))
          anc.state <- ancState(phy, y)
          
          colnames(output.mat) <- 1:dim(output.mat)[2]
          colnames(output.mat)[1:5] <- c("lnL", "AIC", "AICc", "sigma.sq", "anc.state")
          colnames(output.mat)[6:(6 + nSplits)] <- paste0("rates", 1:(nSplits + 1))
          start.time <- 7 + nSplits
          colnames(output.mat)[start.time:(start.time + nSplits - 1)] <- paste0("time.split", 1:nSplits)
          output.mat[1, 1:5] <- c(log.lik, aic, aicc, as.numeric(bm.model[[1]]), anc.state)
          print("BM model", quote = F)
          print(output.mat[1, 1:5])
          fixed.time <- NULL
          vo.out <- list()
          
          for (q in 1:nSplits) {
            all.models <- sapply(splitTime , function(x_times) {
              splitTimeInt <- sort(c(x_times, fixed.time))
              nSplitInt.n <- length(splitTimeInt) + 1
              rateVec <- rep(1, nSplitInt.n)
              var.timeslice <- function(rate.vec) transformPhylo.ll(y, phy, timeRates = rate.vec, splitTime = splitTimeInt, model = "timeSlice", meserr = meserr, covPIC = covPIC)[[2]]
              vo <- optim(rateVec, var.timeslice, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
              c(vo$par, vo$value)
            })
            best.model.n <- which.max(all.models[3, ])
            shift.time <- splitTime[best.model.n]
            fixed.time <- c(fixed.time, shift.time)
            splitTime <- splitTime[-which(splitTime == shift.time)]
            log.lik <- tail(all.models[, best.model.n], 1)
            rates <- all.models[, best.model.n][-dim(all.models)[1]]
            phy.temp <- transformPhylo(y = y, phy = phy, timeRates = rates , splitTime = fixed.time, model = "timeSlice")
            anc.state <- ancState(phy.temp, y)
            bVar <- as.numeric(likTraitPhylo(y, phy.temp, covPIC = covPIC)[[1]])
            param <- length(rates) + 2
            aic <- aic.fun(log.lik, param)
            aicc <- aicc.fun(log.lik, param, Ntip(phy))
            output.mat[q + 1, 1:5] <- c(log.lik, aic, aicc, bVar, anc.state)
            output.mat[q + 1, 6:(5 + length(rates))] <- rates
            start.time <- tail(1:dim(output.mat)[2], nSplits)[1:q]
            output.mat[q + 1, start.time] <- sort(fixed.time)
            na.test <- which(is.na(output.mat[q + 1, ]))
            print(paste("shift", q), quote = F)
            na.test <- which(is.na(output.mat[q + 1, ]))
            if (length(na.test) > 0) {
              print(output.mat[q + 1,-which(is.na(output.mat[q + 1, ]))])
            } else {
              print(output.mat[q + 1,])
            }
            rateVec <- rep(1, length(rates))
            var.timeslice <- function(rate.vec)
                transformPhylo.ll(y, phy, timeRates = rate.vec, splitTime = fixed.time, model = "timeSlice", meserr = meserr, covPIC = covPIC)[[2]]
            vo.out[[q]] <- optim(rateVec, var.timeslice, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
          }
        } else {
          nSplits <- length(splitTime) + 1
          rateVec <- rep(1, nSplits)
          output.mat <- matrix(NA, nrow = 2, ncol = 6 +  2 * length(splitTime))
          colnames(output.mat) <- 1:dim(output.mat)[2]
          colnames(output.mat)[1:5] <- c("lnL", "AIC", "AICc", "sigma.sq", "anc.state")
          colnames(output.mat)[6:(5 + length(splitTime)  + 1)] <- paste0("rate", 1:(length(splitTime)  + 1))
          start.time <- 7 + length(splitTime)
          colnames(output.mat)[start.time:(start.time + nSplits - 2)] <- paste0("time.split", 1:(length(splitTime)))
          bm.model <- likTraitPhylo(y, phy, covPIC = covPIC)
          log.lik <- as.numeric(bm.model[[2]])
          aic <- aic.fun(log.lik, 2)
          aicc <- aicc.fun(log.lik, 2, Ntip(phy))
          anc.state <- ancState(phy, y)
          output.mat[1, 1:5] <- c(log.lik, aic, aicc, as.numeric(bm.model[[1]]), anc.state)
          print("BM model")
          print(output.mat[1, 1:5])
          var.timeslice <- function(rate.vec) transformPhylo.ll(y, phy, timeRates = rate.vec, splitTime = splitTime, model = "timeSlice")[[2]]
          vo.out <- optim(rateVec, var.timeslice, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)
          log.lik <- vo.out$value
          rates <- vo.out$par
          param <- length(rates) + 2
          aic <- aic.fun(log.lik, param)
          aicc <- aicc.fun(log.lik, param, Ntip(phy))
          phyTemp <- transformPhylo(y = y, phy = phy, timeRates = rates, splitTime = splitTime, model = "timeSlice")
          bVar <- as.numeric(likTraitPhylo(y, phyTemp, covPIC = covPIC)[[1]])
          anc.state <- ancState(phyTemp, y)
          output.mat[2, 1:5] <- c(log.lik, aic, aicc, bVar, anc.state)
          output.mat[2, 6:(5 + length(rates))] <- rates
          start.rates <- 6 + length(rates)
          start.time <- tail(1:dim(output.mat)[2], nSplits - 1)
          output.mat[2, start.time] <- splitTime
          print("shiftModel")
          print(output.mat[2,])
        }
       	out <- list()
        out$timeSlice <- output.mat
        out$optim.output <- vo.out
        out$phy <- phy
        out$y <- y
      }
    )
    return(out)
  }