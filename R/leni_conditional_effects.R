#' @name leni_conditional_effects
#' @keywords internal
#' @noRd

leni_conditional_effects <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = FALSE,
    bootSeed = NULL,
    model.class = "lme",
    data = NULL,
    coef.idx = NULL,
    modx = NULL,
    modx.idx = NULL,
    ci = ci,
    ...
){

  # Select Relevant Transformations
  expr <- if(target_fx == "quadratic"){
    list(
      f_a0 = expression(b0),
      f_a0_modx = expression(b3),
      f_ax = expression(-b1/(2 * b2)),
      f_ax_modx = expression((-(b1 + b4)/(2 * (b2 + b5)))-(-b1/(2 * b2))),
      f_ay = expression(b0 - (b1^2/(4 * b2))),
      f_ay_modx = expression(((b0 + b3) - (b1 + b4)^2/(4 * (b2 + b5))) -
                               (b0 - (b1^2/(4 * b2)))),
      f_g = expression(-(b1^2/(4 * b2))),
      f_g_modx = expression((-(b1 + b4)^2/(4 * (b2 + b5))) -
                              (-(b1^2/(4 * b2)))),
      f_ac = expression(b2),
      f_ac_modx = expression(b5)
    )
  } else if(target_fx == "cubic"){
    list(
      f_xn = expression(-b2 / (3 * b3)),
      f_xn_modx = expression((-(b2 + b6)/(3 * (b3 + b7))) - (-b2 / (3 * b3))),
      f_yn = expression(b0 - (b1 * b2) / (3 * b3) + (2 * b2^3) / (27 * b3^2)),
      f_yn_modx = expression(
        ((b0 + b4) - ((b1 + b5) + (b2 + b6))/(3 * (b3 + b7)) + (2 * (b2 + b6)^3)/(27 * (b3 + b7)^2)) -
          (b0 - (b1 * b2) / (3 * b3) + (2 * b2^3) / (27 * b3^2))
      ),
      f_d  = expression(sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))),
      f_d_modx = expression(
        sqrt(((b2 + b6)^2 - 3 * (b3 + b7) * (b1 + b5)) / (9 * (b3 + b7)^2)) -
          sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))
      ),
      f_h  = expression(-2 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^3),
      f_h_modx = expression(
        (-2 * (b3 + b7) * sqrt(((b2 + b6)^2 - 3 * (b3 + b7) * (b1 + b5)) / (9 * (b3 + b7)^2))^3) -
          (-2 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^3)
      ),
      f_bn = expression(-3 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^2),
      f_bn_modx = expression(
        (-3 * (b3 + b7) * sqrt(((b2 + b6)^2 - 3 * (b3 + b7) * (b1 + b5)) / (9 * (b3 + b7)^2))^2) -
          (-3 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^2)
      )
    )
  }
  expr <- expr[c(paste0("f_", tolower(theta)),
                 paste0("f_", tolower(theta),"_modx"))]

  # Save Out Relevant Parameter Estimates
  if(grepl(model.class,"lm")){
    crit.val <- stats::qt((1 - ci) / 2,
                          df = nrow(model$model),
                          lower.tail = FALSE)
  } else if(grepl(model.class,"lme")){
    crit.val <- stats::qt((1 - ci) / 2,
                   df = length(unique(model@flist[[1]])),
                   lower.tail = FALSE)
  }

  # Run Bootstrap or Analytic Transformations
  if(is.numeric(bootstrap) & bootstrap > 0){
    if(!is.null(bootSeed)){set.seed(bootSeed)}
    bootResults <- boot::boot(
      data = if (grepl(model.class,"lm")){
        data
      } else {
        data |> dplyr::group_nest(model@flist[[1]])
      },
      statistic = conditional_effects_bootstrap,
      R = bootstrap,
      model = model,
      model.class = model.class,
      target_fx = target_fx,
      expr = expr,
      coef.idx = coef.idx,
      modx.idx = modx.idx,
      ...)

    return(
      list(
        fixed_effects = data.frame(
          row.names = c(theta, paste0(theta, "_", modx)),
          est = bootResults$t0,
          se = apply(bootResults$t, 2, function(x) stats::sd(x, na.rm = TRUE)),
          ci.lower = tryCatch(
            expr = {sapply(1:ncol(bootResults$t),function(x){boot::boot.ci(bootResults, type = "bca", index = x)$bca[4]})},
            error = function(e){bootResults$t0 - crit.val*apply(bootResults$t, 2, function(x) stats::sd(x, na.rm = TRUE))}
          ),
          ci.upper = tryCatch(
            expr = {sapply(1:ncol(bootResults$t),function(x){boot::boot.ci(bootResults, type = "bca", index = x)$bca[5]})},
            error = function(e){bootResults$t0 + crit.val*apply(bootResults$t, 2, function(x) stats::sd(x, na.rm = TRUE))}
          ),
          est.robust = apply(bootResults$t, 2, function(x){
            if (stats::sd(x, na.rm = TRUE) > 0) bayestestR::map_estimate(stats::na.omit(x))
            else mean(x, na.rm = TRUE)
          }),
          se.robust = apply(bootResults$t, 2, function(x) trim_sd(x, na.rm = TRUE))
        ),
        bootstrap_samples = bootResults
      )
    )
  } else {
    if (grepl(target_fx, "quadratic")){
      c(b0,b1,b2) %tin% read.coefs(model)[coef.idx]
      c(b3,b4,b5) %tin% read.coefs(model)[modx.idx]
    } else if (grepl(target_fx, "cubic")){
      c(b0,b1,b2,b3) %tin% read.coefs(model)[coef.idx]
      c(b4,b5,b6,b7) %tin% read.coefs(model)[modx.idx]
    }
    nl_theta <- sapply(expr, eval, envir = environment())

    J <- matrix(c(
      sapply(expr, function(x){
        sapply(paste0("b",0:(length(c(coef.idx,modx.idx))-1)),
               function(y){eval(stats::D(x,y))})})
    ),
    nrow = length(c(coef.idx,modx.idx)),
    ncol = length(c(coef.idx,modx.idx)),
    byrow = FALSE)

    return(
      list(
        fixed_effects = data.frame(
          row.names = names(nl_theta),
          est = nl_theta,
          se = sqrt(diag(t(J) %*% stats::vcov(model)[c(coef.idx,modx.idx),
                                              c(coef.idx,modx.idx)] %*% J)),
          ci.lower = nl_theta -
            crit.val*sqrt(diag(t(J) %*% stats::vcov(model)[c(coef.idx,modx.idx),
                                                    c(coef.idx,modx.idx)] %*% J)),
          ci.upper = nl_theta +
            crit.val*sqrt(diag(t(J) %*% stats::vcov(model)[c(coef.idx,modx.idx),
                                                    c(coef.idx,modx.idx)] %*% J))
        ),
        ACOV_theta = structure(
          as.matrix(t(J) %*% stats::vcov(model)[c(coef.idx,modx.idx),
                                         c(coef.idx,modx.idx)] %*% J),
          dimnames = list(substr(names(nl_theta), 3, nchar(names(nl_theta))),
                          substr(names(nl_theta), 3, nchar(names(nl_theta))))
        )
      )
    )
  }
}

conditional_effects_bootstrap <- function(data, indices, model,
                                          model.class, target_fx, expr,
                                          coef.idx, modx.idx){
  if (grepl(model.class,"lm")){
    fit <- stats::update(model, data = data[indices,])
  } else if (grepl(model.class,"lme")){
    fit <- stats::update(model, data = data[indices,] |> tidyr::unnest(cols=c(data)))
  }

  if (grepl(target_fx, "quadratic")){
    c(b0,b1,b2) %tin% read.coefs(fit)[coef.idx]
    c(b3,b4,b5) %tin% read.coefs(fit)[modx.idx]
  } else if (grepl(target_fx, "cubic")){
    c(b0,b1,b2,b3) %tin% read.coefs(fit)[coef.idx]
    c(b4,b5,b6,b7) %tin% read.coefs(fit)[modx.idx]
  }

  nl_theta <- sapply(expr, eval, envir = environment())
  return(
    stats::setNames(nl_theta, substr(names(nl_theta), 3, nchar(names(nl_theta))))
  )
}
