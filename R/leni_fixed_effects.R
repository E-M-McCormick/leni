leni_fixed_effects <- function(
    model = NULL,
    model.class = NULL,
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
  expr <- expr[paste0("f_", theta)]

  # Save Out Relevant Parameter Estimates
  read.coefs <- if(model.class == "lm"){
    stats::coef
  } else if(model.class == "lme"){
    lme4::fixef
  }

  # Run bootstrap
  if(is.numeric(bootstrap) & bootstrap > 0){
    results <- boot::boot(
      data = if(model.class == "lm"){dat
             } else {dat |> dplyr::group_nest(names(model@flist))},
      statistic = if(model.class == "lm"){lm_fixed_effects_bootstrap
                  } else {lme_fixed_effects_bootstrap},
      R = bootstrap,
      model.obj = model)
  } else {
    invisible(
      mapply(assign,
             paste0("b",0:(length(read.coefs(model))-1)),
             read.coefs(model),
             MoreArgs = list(envir = parent.frame()))
    )
    nl_theta <- sapply(expr, eval)
    point_estimates <- setNames(nl_theta,
                                substr(names(nl_theta),
                                       3,
                                       nchar(names(nl_theta))))
    J <- matrix(c(
      sapply(expr, function(x){
        sapply(paste0("b",0:(length(read.coefs(model))-1)),
               function(y){eval(D(x,y))})})
      ), nrow = length(read.coefs(model)), ncol = length(read.coefs(model)),
      byrow = FALSE)
    ACOV_theta <- t(J) %*% vcov(model) %*% J

  }
}

bootstrapCall <- function(call = NULL){
  call_str <- call |> deparse1(width.cutoff = 500, collapse = "")
  if(grepl("data.*=.*?,", call_str)){
    return(
      call_str |>
        gsub("data.*=.*?,",
             "data = dat[indices,] |> tidyr::unnest(cols=c(data)),",
             x = _) |>
        str2lang()
    )
  } else {
    return(
      call_str |>
        sub(",.*?,",
            ", data = dat[indices,] |> tidyr::unnest(cols=c(data)),",
            x = _) |>
        str2lang()

    )
  }
} # need to figure out something about getting package calls (e.g., lmerTest::lmer) (maybe another function to assign main function out of the relevant package)
# add in functionality for glmmTMB package

lm_fixed_effects_bootstrap <- function(data, indices, model.obj){
  fit <- eval(boostrapCall(getCall(model.obj)))
  invisible(
    mapply(assign,
         paste0("b",0:(length(read.coefs(fit))-1)),
         read.coefs(fit),
         MoreArgs = list(envir = parent.frame()))
  )
  nl_theta <- sapply(expr, eval)
  return(
    setNames(nl_theta, substr(names(nl_theta), 3, nchar(names(nl_theta))))
  )
}
lme_fixed_effects_bootstrap <- function(data, indices, model.obj){
  fit <- eval(boostrapCall(getCall(model.obj)))
  invisible(
    mapply(assign,
           paste0("b",0:(length(read.coefs(fit))-1)),
           read.coefs(fit),
           MoreArgs = list(envir = parent.frame()))
  )
  nl_theta <- sapply(expr, eval)
  return(
    setNames(nl_theta, substr(names(nl_theta), 3, nchar(names(nl_theta))))
  )
}
