# Optimisation Criterion: pooled skewness as defined in Rojas-Perillas et al. (2020)

#' Pooled skewness of the area-specific random effect and error term
#' 
#' Sollte man sich auch noch die Wölbung ansehen? Jede um den MW symmetrische
#' Verteilung besitzt Schiefe 0, nicht nur die NV (s.a. MH Th. 3.20).
#'
#' @param fixed two-sided linear formula object which is passed to `nlme::lme()`
#' @param smp_data Sample data
#' @param smp_domains Sample domains
#' @param transformation The data-driven transformaton applied to the sample data
#' @param lambda The transformation parameter rela5tinv to the data-driven transformation
#'
#' @return The pooled skewness as given in Rojas-Perillas et al. (2020); double.
#' @export
#'
#' @examples

pooled_skewness <- function(fixed, smp_data, smp_domains, transformation, lambda) {
	# Applying the transformation to the sample data 
	# vmtl am besten mit emdi:::std_data_transformed
	# Estimating the linear model to obtain the required variance components
	smp_data[as.character(fixed[[2]])] <- emdi:::box_cox(
		y = smp_data[as.character(fixed[[2]])]
		, lambda = lambda
		, shift = 0)$y
	
	model_REML <- NULL
	try(
		model_REML <- nlme::lme(fixed = fixed
											, data = smp_data
											, random = as.formula(paste0("~ 1 | as.factor(", smp_domains, ")"))
											, method = "REML"
											, keep.data = FALSE
		)
		, silent = TRUE
	)
	
	if (is.null(model_REML)) {
		stop(strwrap(prefix = " ", initial = ""
								 , "The likelihood does not converge. One reason could be that
								 the interval for the estimation of an optimal transformation
								 parameter is not appropriate. Try another interval. See also
								 help(ebp)."))
	}
	else {
		model_REML <- model_REML
	}
	# Extract residuals
	resid_u <- resid(model_REML, type = "response", level = 0) # Level 1/area level 
	resid_e <- resid(model_REML, type = "response", level = 1) # Level 1/individual level\
	
	# Extract variance components
	sigma2_u_est <- as.numeric(nlme::VarCorr(model_REML)[1, 1])
	sigma2_e_est <- as.numeric(nlme::VarCorr(model_REML)[2, 1])
	w <- sigma2_e_est / (sigma2_e_est + sigma2_u_est)
	
	# Compute the empirical skewness The empirical skewness is computed/estimated
	# without the Bessel correction. R's var()/sd() use n-1 in the denominator.
	skewness <- function(x) {
		(sum((x - mean(x))^3) / length(x)) / (sd(x) * (length(x) - 1) / length(x))^3
	}
	
	skew_u <- skewness(resid_u)
	skew_e <- skewness(resid_e)
	
	# Pooled skewness
	w * abs(skew_e) + (1 - w) * abs(skew_u)
}
