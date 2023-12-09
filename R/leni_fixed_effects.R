leni_fixed_effects <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = FALSE,
    model.class = "lme",
    data = NULL,
    coef.idx = NULL,
    ...
){

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

  # Save Out Relevant Parameter Estimates
  read.coefs <- if(model.class == "lm"){
    stats::coef
  } else if(model.class == "lme"){
    lme4::fixef
  }

  # Run bootstrap
  if(is.numeric(bootstrap) & bootstrap > 0){
    if(!is.null(bootSeed)){set.seed(bootSeed)}
    bootResults <- boot::boot(
      data = if (model.class == "lm"){
               data
             } else {
               data |> dplyr::group_nest(names(model@flist))
             },
      statistic = if (model.class == "lm"){
                    lm_fixed_effects_bootstrap
                  } else {
                    lme_fixed_effects_bootstrap
                  },
      R = bootstrap,
      model = model,
      ...)

    return(
      list(
        fixed_effects = bootResults$t0,
        fixed_effects_se = apply(bootResults$t,
                                 2,
                                 function(x) sd(x, na.rm = TRUE)),
        fixed_effects_robust = apply(bootResults$t,
                                     2,
                                     function(x)
                                       bayestestR::map_estimate(
                                         na.omit(x))
                                     ),
        fixed_effects_se_robust = apply(bootResults$t,
                                        2,
                                        function(x) IQR(x, na.rm = TRUE)),
        bootstrap_samples = bootResults
      )
    )
  } else {
    if (grepl(target_fx, "quadratic")){
      c(b0,b1,b2) %tin% read.coefs(model)[1:3]
    } else if (grepl(target_fx, "cubic")){
      c(b0,b1,b2,b3) %tin% read.coefs(model)[1:4]
    }
    nl_theta <- sapply(expr, eval, envir = environment())

    J <- matrix(c(
      sapply(expr, function(x){
        sapply(paste0("b",0:(length(read.coefs(model))-1)),
               function(y){eval(D(x,y))})})
      ), nrow = length(read.coefs(model)), ncol = length(read.coefs(model)),
      byrow = FALSE)

    return(
      list(
        fixed_effects = setNames(nl_theta,
                                 substr(names(nl_theta),
                                        3,
                                        nchar(names(nl_theta)))),
        fixed_effects_se = sqrt(diag(t(J) %*% vcov(model) %*% J)),
        ACOV_theta = t(J) %*% vcov(model) %*% J
      )
    )
  }
}


lm_fixed_effects_bootstrap <- function(data, indices, model, ...){
  d <- data[indices,]
  fit <- update(model, data = d)

  c(b0,b1,b2) %tin% read.coefs(fit)[1:3]
  if(grepl(target_fx, "cubic")){b3 <- read.coefs(fit)[4]}

  nl_theta <- sapply(expr, eval, envir = environment())
  return(
    setNames(nl_theta, substr(names(nl_theta), 3, nchar(names(nl_theta))))
  )
}

lme_fixed_effects_bootstrap <- function(data, indices, model, ...){
  d <- data[indices,]
  fit <- update(model, data = d |> tidyr::unnest(cols=c(data)))

  c(b0,b1,b2) %tin% read.coefs(fit)[1:3]
  if(grepl(target_fx, "cubic")){b3 <- read.coefs(fit)[4]}

  nl_theta <- sapply(expr, eval, envir = environment())
  return(
    setNames(nl_theta, substr(names(nl_theta), 3, nchar(names(nl_theta))))
  )
}
