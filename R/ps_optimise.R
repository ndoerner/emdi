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
  weighted_residuals <- function(mp = point_estim$model_par, f = framework, fixed = fixed, lambda = lambda) {
    split(x = framework$smp_data, f = framework$smp_domains_vec) |>
      lapply(function(dom) {
        d <- which(unique(f$smp_domains_vec) %in% unique(dom[[f$smp_domains]]))
        transformed_y <- box_cox_std(dom[, all.vars(fixed)[1]], lambda)
        covariates <- dom[all.vars(fixed)[-1]]
        factors <- vapply(covariates, is.factor, FUN.VALUE = logical(1))
        covariates[factors] <- vapply(covariates[factors], function(x) {
          as.numeric(x) - 1
        }, FUN.VALUE = numeric(nrow(covariates)))
        covariates <- as.matrix(cbind(1, covariates))

        # Summe aus der Definition nicht relevant, da Skalarprodukt/Punktprodukt
        x_iw <- sum(dom$weight)^-1 * dom$weight %*% covariates
        y_iw <- sum(dom$weight)^-1 * dom$weight %*% transformed_y

        fitted <- as.vector(covariates %*% mp$betas) + as.vector(mp$gammaw[d] * (y_iw  - x_iw %*% mp$betas))

        # Return Pearson residuals. sigma_w_2 is based on G. (16)
        sigma_w_2 <- mp$sigmau2est * (1 - mp$gammaw[d]) + mp$sigmae2est
        (transformed_y - fitted) / sqrt(sigma_w_2)
      }) |> unsplit(f = framework$smp_domains_vec)
  }
  resid <- weighted_residuals(point_estim$model_par, framework, fixed, lambda) # resid(point_estim$model, level = 0, type = "pearson")
  ranef <- point_estim$model_par$rand_eff[framework$dist_obs_dom]

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