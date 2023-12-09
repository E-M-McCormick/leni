leni_random_effects <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = 5000,
    model.class = "lme",
    data = NULL,
    coef.idx = NULL,
    ...){

  # Select Relevant Transformations
  expr <- if(target_fx == "quadratic"){
    list(
      f_a0 = expression(b0),
      f_ax = expression(-b1/(2 * b2)),
      f_ay = expression(b0 - (b1^2/(4 * b2))),
      f_g = expression(- (b1^2/(4 * b2))),
      f_ac = expression(b2)
    )
  } else if(target_fx == "cubic"){
    list(
      f_xn = expression(-b2 / (3 * b3)),
      f_yn = expression(b0 - (b1 * b2) / (3 * b3) + (2 * b2^3) / (27 * b3^2)),
      f_d  = expression(sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))),
      f_h  = expression(-2 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^3),
      f_bn = expression(-3 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^2)
    )
  }
  expr <- expr[paste0("f_", tolower(theta))]

  if (grepl(target_fx, "quadratic")){
    c(b0,b1,b2) %tin% lme4::fixef(model)[1:3]
  } else if (grepl(target_fx, "cubic")){
    c(b0,b1,b2,b3) %tin% lme4::fixef(model)[1:4]
  }

  J <- matrix(c(
    sapply(expr, function(x){
      sapply(paste0("b",0:(length(lme4::fixef(model))-1)),
             function(y){eval(D(x,y))})})
  ), nrow = length(lme4::fixef(model)), ncol = length(lme4::fixef(model)),
  byrow = FALSE)

  # Save Out Relevant Parameter Estimates
  vcor_matrix <- matrix(0, nrow = length(theta), ncol = length(theta))
  vcor_matrix[1:nrow(lme4::VarCorr(model)$id),
              1:ncol(lme4::VarCorr(model)$id)] <- lme4::VarCorr(model)$id
  tau_theta <- t(J) %*% vcor_matrix %*% J

  if(is.numeric(bootstrap) & bootstrap > 0){
    if(!is.null(bootSeed)){set.seed(bootSeed)}
    bootResults <- boot::boot(
      data = data |> dplyr::group_nest(names(model@flist)),
      statistic = lme_random_effects_bootstrap,
      R = bootstrap,
      model = model,
      ...)

    return(
      list(
        tau_theta = tau_theta,
        tau_theta_bootstrap = apply(
          bootResults$t,
          2,
          function(x) mean(x, na.rm = TRUE)),
        tau_theta_bootstrap_se = apply(
          bootResults$t,
          2,
          function(x) sd(x, na.rm = TRUE)),
        tau_theta_bootstrap_robust = apply(
          bootResults$t,
          2,
          function(x){
            if (sd(x, na.rm = TRUE) > 0){
              bayestestR::map_estimate(na.omit(x))
            } else {
              0
            }
          }),
        tau_theta_bootstrap_se_robust = apply(
          bootResults$t,
          2,
          function(x) IQR(x, na.rm = TRUE)),
        bootstrap_samples = bootResults
      )
    )
  }
}

lme_random_effects_bootstrap <- function(data, indices, model, ...){
  d <- data[indices,]
  fit <- update(model, data = d |> tidyr::unnest(cols=c(data)))

  c(b0,b1,b2) %tin% lme4::fixef(fit)[1:3]
  if(grepl(target_fx, "cubic")){b3 <- lme4::fixef(fit)[4]}
  J <- matrix(c(
    sapply(expr, function(x){
      sapply(paste0("b",0:(length(lme4::fixef(fit))-1)),
             function(y){eval(D(x,y))})})
  ), nrow = length(lme4::fixef(fit)), ncol = length(lme4::fixef(fit)),
  byrow = FALSE)

  vcor_matrix <- matrix(0, nrow = length(theta), ncol = length(theta))
  vcor_matrix[1:nrow(lme4::VarCorr(fit)$id),
              1:ncol(lme4::VarCorr(fit)$id)] <- lme4::VarCorr(fit)$id
  tau_theta <- t(J) %*% vcor_matrix %*% J

  nl_theta <- sapply(expr, eval, envir = environment())
  return(
    tau_theta[lower.tri(tau_theta, diag = TRUE)]
  )
}
