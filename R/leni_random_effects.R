leni_random_effects <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = FALSE,
    bootSeed = NULL,
    model.class = "lme",
    data = NULL,
    coef.idx = NULL,
    ci = ci,
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
    crit.val <- NULL                              # FIX ME
  } else if (model.class == "lme"){
    crit.val <- qt((1 - ci) / 2,
                   df = length(unique(leni_realdata@flist[[1]])),
                   lower.tail = FALSE)
  }

  if (grepl(target_fx, "quadratic")){
    c(b0,b1,b2) %tin% read.coefs(model)[coef.idx]
    cov_theta_names <- c(
      theta[1], paste0(theta[1],theta[2]), paste0(theta[1],theta[3]),
      theta[2], paste0(theta[2],theta[3]), theta[3])
  } else if (grepl(target_fx, "cubic")){
    c(b0,b1,b2,b3) %tin% read.coefs(model)[coef.idx]
    cov_theta_names <- c(
      theta[1], paste0(theta[1],theta[2]), paste0(theta[1],theta[3]),
      paste0(theta[1],theta[4]), theta[2], paste0(theta[2],theta[3]),
      paste0(theta[2],theta[4]), theta[3], paste0(theta[3],theta[4]), theta[4])
  }

  J <- matrix(c(
    sapply(expr, function(x){
      sapply(paste0("b",0:(length(coef.idx)-1)),
             function(y){eval(D(x,y))})})
  ), nrow = length(coef.idx),
  ncol = length(coef.idx),
  byrow = FALSE)

  #  Run Analytic Transformations
  vcor_matrix <- matrix(0, nrow = length(coef.idx), ncol = length(coef.idx))
  vcor_matrix[1:nrow(lme4::VarCorr(model)$id),
              1:ncol(lme4::VarCorr(model)$id)] <- lme4::VarCorr(model)[[1]]

  # Run Bootstrap Transformations
  if(is.numeric(bootstrap) & bootstrap > 0){
    if(!is.null(bootSeed)){set.seed(bootSeed)}
    bootResults <- boot::boot(
      data = data |> dplyr::group_nest(model@flist[[1]]),
      statistic = random_effects_bootstrap,
      R = bootstrap,
      model = model,
      model.class = model.class,
      target_fx = target_fx,
      expr = expr,
      coef.idx = coef.idx,
      ...)

    return(
      list(
        tau_theta = structure(
          as.matrix(t(J) %*% vcor_matrix %*% J),
          dimnames = list(theta, theta)
        ),
        random_effects = data.frame(
          row.names = cov_theta_names,
          boot.est = apply(bootResults$t, 2, function(x) mean(x, na.rm = TRUE)),
          boot.se = apply(bootResults$t, 2, function(x) sd(x, na.rm = TRUE)),
          boot.ci.lower = tryCatch(
            expr = {sapply(1:ncol(bootResults$t), function(x){boot::boot.ci(bootResults, type = "bca", index = x)$bca[4]})},
            error = function(e){apply(bootResults$t, 2, function(x) mean(x, na.rm = TRUE)) -
                crit.val*apply(bootResults$t, 2, function(x) sd(x, na.rm = TRUE))}
          ),
          boot.ci.upper = tryCatch(
            expr = {sapply(1:ncol(bootResults$t),function(x){boot::boot.ci(bootResults, type = "bca", index = x)$bca[5]})},
            error = function(e){
              apply(bootResults$t, 2, function(x) mean(x, na.rm = TRUE)) +
                crit.val*apply(bootResults$t, 2, function(x) sd(x, na.rm = TRUE))}
          ),
          boot.est.robust = apply(bootResults$t, 2, function(x){
            if (sd(x, na.rm = TRUE) > 0) bayestestR::map_estimate(na.omit(x))
            else mean(x, na.rm = TRUE)
          }),
          boot.se.robust = apply(bootResults$t, 2, function(x) IQR(x, na.rm = TRUE))
        ),
        bootstrap_samples = bootResults
      )
    )
  } else {
    return(list(
      tau_theta = structure(
        as.matrix(t(J) %*% vcor_matrix %*% J),
        dimnames = list(theta, theta))))
  }
}

random_effects_bootstrap <- function(data, indices, model,
                                     model.class, target_fx, expr, coef.idx){
  fit <- update(model, data = data[indices,] |> tidyr::unnest(cols=c(data)))

  c(b0,b1,b2) %tin% lme4::fixef(fit)[1:3]
  if(grepl(target_fx, "cubic")){b3 <- lme4::fixef(fit)[4]}
  J <- matrix(c(
    sapply(expr, function(x){
      sapply(paste0("b",0:(length(coef.idx)-1)),
             function(y){eval(D(x,y))})})
  ),
  nrow = length(coef.idx),
  ncol = length(coef.idx),
  byrow = FALSE)

  vcor_matrix <- matrix(0, nrow = length(coef.idx), ncol = length(coef.idx))
  vcor_matrix[1:nrow(lme4::VarCorr(fit)$id),
              1:ncol(lme4::VarCorr(fit)$id)] <- lme4::VarCorr(fit)$id
  tau_theta <- t(J) %*% vcor_matrix %*% J

  nl_theta <- sapply(expr, eval, envir = environment())
  return(
    tau_theta[lower.tri(tau_theta, diag = TRUE)]
  )
}
