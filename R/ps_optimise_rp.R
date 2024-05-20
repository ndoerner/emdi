ps_optimise_rp <- function(lambda
                           , fixed
                           , transformation = "box.cox"
                           , framework
                           , L
                           , keep_data = FALSE
) {
  dat <- framework$smp_data
  y <- all.vars(fixed)[1]
  dat[y] <- box_cox(y = dat[y], lambda = lambda)$y
  mod <- NULL
  try(
    mod <- lme(fixed = fixed # as.formula(paste0(fixed, "- 1"))
               , data = dat
               , random = as.formula(paste0("~ 1 | as.factor(", framework$smp_domains, ")"))
               , method = "REML"
    )
  )
  if (is.null(mod)) return(999999)

  resid <- resid(mod, level = 0, type = "pearson")
  ranef <- as.matrix(random.effects(mod))[, 1] # point_estim$model_par$rand_eff[framework$dist_obs_dom]

  sigma2_e_est <- mod$sigma^2
  sigma2_u_est <- as.numeric(VarCorr(mod)[1, 1])
  w <- sigma2_e_est / (sigma2_e_est + sigma2_u_est)

  # Compute the empirical skewness The empirical skewness is computed/estimated
  # without the Bessel correction. R's var()/sd() use n-1 in the denominator.
  skewness <- function(x) {
    (sum((x - mean(x))^3)/ length(x)) / (sum((x - mean(x))^2) / length(x))^(3/2)
  }

  skew_e <- skewness(resid)
  skew_u <- skewness(ranef)

  # Pooled skewness
  w * abs(skew_e) + (1 - w) * abs(skew_u)
}