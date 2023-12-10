#' @name leni_fixed_effects
#' @keywords internal
#' @noRd

leni_fixed_effects <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = FALSE,
    bootSeed = NULL,
    model.class = "lme",
    data = NULL,
    coef.idx = NULL,
    ci = 0.95,
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
  if (model.class == "lm"){
    crit.val <- nrow(model$model)
  } else if (model.class == "lme"){
    crit.val <- stats::qt((1 - ci) / 2,
                   df = length(unique(model@flist[[1]])),
                   lower.tail = FALSE)
  }

  # Run Bootstrap or Analytic Transformations
  if(is.numeric(bootstrap) & bootstrap > 0){
    if(!is.null(bootSeed)){set.seed(bootSeed)}
    bootResults <- boot::boot(
      data = if (model.class == "lm"){
          data
        } else {
          data |> dplyr::group_nest(model@flist[[1]])
        },
      statistic = fixed_effects_bootstrap,
      R = bootstrap,
      model = model,
      model.class = model.class,
      target_fx = target_fx,
      expr = expr,
      coef.idx = coef.idx,
      ...)

    return(
      list(
        fixed_effects = data.frame(
          row.names = theta,
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
          se.robust = apply(bootResults$t, 2, function(x) stats::IQR(x, na.rm = TRUE))
        ),
        bootstrap_samples = bootResults
      )
    )
  } else {
    if (grepl(target_fx, "quadratic")){
      c(b0,b1,b2) %tin% read.coefs(model)[coef.idx]
    } else if (grepl(target_fx, "cubic")){
      c(b0,b1,b2,b3) %tin% read.coefs(model)[coef.idx]
    }
    nl_theta <- sapply(expr, eval, envir = environment())

    J <- matrix(c(
      sapply(expr, function(x){
        sapply(paste0("b",0:(length(read.coefs(model)[coef.idx])-1)),
               function(y){eval(stats::D(x,y))})})
      ),
      nrow = length(read.coefs(model)[coef.idx]),
      ncol = length(read.coefs(model)[coef.idx]),
      byrow = FALSE)

    return(
      list(
        fixed_effects = data.frame(
          row.names = theta,
          est = nl_theta,
          se = sqrt(diag(t(J) %*% stats::vcov(model)[coef.idx,coef.idx] %*% J)),
          ci.lower = nl_theta -
            crit.val*sqrt(diag(t(J) %*% stats::vcov(model)[coef.idx,coef.idx] %*% J)),
          ci.upper = nl_theta +
            crit.val*sqrt(diag(t(J) %*% stats::vcov(model)[coef.idx,coef.idx] %*% J))
        ),
        ACOV_theta = structure(
          as.matrix(t(J) %*% stats::vcov(model)[coef.idx,coef.idx] %*% J),
          dimnames = list(substr(names(nl_theta), 3, nchar(names(nl_theta))),
                          substr(names(nl_theta), 3, nchar(names(nl_theta))))
        )
      )
    )
  }
}


fixed_effects_bootstrap <- function(data, indices, model,
                                    model.class, target_fx, expr, coef.idx){
  if (model.class == "lm"){
    fit <- stats::update(model, data = data[indices,])
  } else if (model.class == "lme"){
    fit <- stats::update(model, data = data[indices,] |> tidyr::unnest(cols = c(data)))
  }

  c(b0,b1,b2) %tin% read.coefs(fit)[coef.idx]
  if(grepl(target_fx, "cubic")){b3 <- read.coefs(fit)[coef.idx[4]]}

  nl_theta <- sapply(expr, eval, envir = environment())
  return(
    stats::setNames(nl_theta, substr(names(nl_theta), 3, nchar(names(nl_theta))))
  )
}
