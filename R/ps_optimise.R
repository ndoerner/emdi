ps_optimise <- function(lambda
                        , fixed
                        , transformation = "box.cox"
                        , framework
                        , L
                        , keep_data = FALSE
) {
  point_estim <- point_estim(
    lambda = lambda,
    framework = framework,
    fixed = fixed,
    transformation = transformation,
    interval = interval,
    L = L,
    keep_data = TRUE
  )

  resid <- resid(point_estim$model, level = 0, type = "pearson")
  ranef <- nlme::ranef(point_estim$model)[[1]]

  sigma2_e_est <- point_estim$model_par$sigmae2est
  sigma2_u_est <- point_estim$model_par$sigmau2est
  w <- sigma2_e_est / (sigma2_e_est + sigma2_u_est)

  # Compute the empirical skewness The empirical skewness is computed/estimated
  # without the Bessel correction. R's var()/sd() use n-1 in the denominator.
  skewness <- function(x) {
    (sum((x - mean(x))^3) / length(x)) / (sd(x) * (length(x) - 1) / length(x))^3
  }

  skew_u <- skewness(resid)
  skew_e <- skewness(ranef)

  # Pooled skewness
  w * abs(skew_e) + (1 - w) * abs(skew_u)
}