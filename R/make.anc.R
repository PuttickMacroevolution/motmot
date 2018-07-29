#' Create design matrix (internal function)
#'
#' This is an internal function to generate the design matrix required to define different means for each hypothesised rate category.
#' @param y The response variable - typically a continuous trait.
#' @param x The explanatory (discrete) variable used to define the hypothesised rate categories. Can be specified as a column number or column name.
#' @param data A data frame containing (minimally) the x and y variables as columns with species names as rownames.
#' @param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates..
#' @return A design matrix
#' @references Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
#' Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas

make.anc <-
function(y, x, data=NULL, common.mean=FALSE) {
	
	if(is.factor(x) == "FALSE"){
		stop("The discrete trait must be a factor")	
	}
	
	m <- model.frame(y ~ x, data)
	x <- model.matrix(y ~ x, m)
	return(x)
}

