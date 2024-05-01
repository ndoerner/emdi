ps_optimise <- function(lambda
                         , fixed
                         , framework
                         , transformation = "box.cox"
                         , L
                         , keep_data = FALSE
) {
  # Data transformation
  transformation_par <- data_transformation(
    fixed = fixed,
    smp_data = framework$smp_data,
    transformation = transformation,
    lambda = lambda
  )
  shift_par <- transformation_par$shift

  # Estimate the mixed model
  mixed_model <- nlme::lme(
    fixed = fixed
    , data = transformation_par$transformed_data
    , random = as.formula(paste0("~ 1 | as.factor(", framework$smp_domains, ")"))
    , method = "REML"
    , keep.data = keep_data
  )
  }