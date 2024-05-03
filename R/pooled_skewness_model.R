#' Pooled skewness of the area-specific random effect and error term
#'
#' @param model A model object describing a linear mixed model
#'
#' @return The pooled skewness as given in Rojas-Perillas et al. (2020); double.
#' @export

pooled_skewness <- function(model) {
  # Extract residuals
  resid_u <- resid(model, type = "response", level = 0) # Level 2/area level;
  resid_e <- resid(model, type = "response", level = 1) # Level 1/individual level\

  # Extract variance components
  sigma2_u_est <- as.numeric(nlme::VarCorr(model)[1, 1])
  sigma2_e_est <- as.numeric(nlme::VarCorr(model)[2, 1])
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