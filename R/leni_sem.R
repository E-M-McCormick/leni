leni_sem <- function(
    target_fx = "cubic",
    theta = c("xN","yN","d","h"),
    start_values = NULL,
    number_obs = 5,
    center_obs = 1,
    spacing = 1,
    verbose = FALSE
){
  # Check NULL arguments
  if(is.null(target_fx)){stop('Target function is not defined. Stopping...')}

  # Define Time Points
  t <- if(length(spacing == 1)){
    seq(1, (spacing*number_obs), by = spacing)
  } else if(length(spacing > 1)) {
    c(1, 1 + cumsum(spacing))
  }

  # Define Target Function & Parameters
  if(grepl(target_fx,"quadratic")){
    tidymessage('Defaulting to the standard alternative quadratic from
                 Cudeck & du Toit, 2002. For the version with gamma
                 please use: target_fx = "quadratic_gamma_a0",
                "quadratic_gamma_ay", or "quadratic_ac"')
    fx <- expression(ay - (ay - a0)*((x_ti / ax) - 1)^2)
    theta <- c("a0","ax","ay")
  } else if(grepl(target_fx,"quadratic_gamma_a0")) {
    fx <- expression(a0 - gamma*(((x_ti / ax) - 1)^2 - 1))
    theta <- c("a0","ax","gamma")
  } else if(grepl(target_fx,"quadratic_gamma_ay")){
    fx <- expression(ay - gamma*((x_ti / ax) - 1)^2)
    theta <- c("ax","ay","gamma")
  } else if(grepl(target_fx,"quadratic_ac")){
    fx <- expression(ay + ac*(x_ti - ax)^2)
    theta <- c("ax","ay","ac")
  } else if(grepl(target_fx,"cubic") |
            grepl(target_fx, "cubic_h")){
    tidymessage('Defaulting to the alternative cubic with h from
                 McCormick, 2024. For the version with betaN
                 please use: target_fx = "cubic_betaN"')
    fx <- expression(yN - (h/2)*(((x_ti-xN)/d)^3 - 3*((x_ti-xN)/d)))
    theta <- c("xN","yN","d","h")
  } else if(grepl(target_fx,"cubic_betaN")){
    fx <- expression(yN - (betaN*d/3)*(((x_ti-xN)/d)^3 - 3*((x_ti-xN)/d)))
    theta <- c("xN","yN","d","betaN")
  } else if(grepl(target_fx, "cubic_multiphase") |
            grepl(target_fx, "cubic_multiphase_h")){
    tidymessage('Defaulting to the alternative multiphase cubic with h from
                 McCormick, 2024. For the version with betaN
                 please use: target_fx = "cubic_multiphase_betaN"')
    fx <- expression(
      yN - (h/2)*(
        (((1/2)*(sqrt((xN-d-x_ti)^2)-sqrt((x_ti-xN-d)^2)))/d)^3 -
          3*(((1/2)*(sqrt((xN-d-x_ti)^2)-sqrt((x_ti-xN-d)^2)))/d))
    )
    theta <- c("xN","yN","d","h")
  } else if(grepl(target_fx, "cubic_multiphase_betaN")){
    fx <- expression(
      yN - (betaN*d/3)*(
        (((1/2)*(sqrt((xN-d-x_ti)^2)-sqrt((x_ti-xN-d)^2)))/d)^3 -
          3*(((1/2)*(sqrt((xN-d-x_ti)^2)-sqrt((x_ti-xN-d)^2)))/d))
    )
    theta <- c("xN","yN","d","betaN")
  } else if(!is.character(target_fx) & !is.expression(target_fx)) {
    stop(tidymessage('Unable to parse target function as a
                     string for conversion to expression.'))
  } else {
    fx <- as.expression(target_fx)
  }
  if(is.null(start_values)){start_values <- rep(0, length(theta))}
  else if(length(start_values) != length(theta)){
    stop(tidymessage('Length of start values does not equal the number of
                     parameters in the target function. Please adjust the input
                     arguments.'))
  }
  names(start_values) <- theta

  # Define Mean Function
  suppressWarnings(mean_fx <- stringi::stri_replace_all_regex(
    fx, theta, paste0(theta,"_m"), vectorize_all = FALSE) |>
    str2expression())

  # Define Item Intercepts and Residuals
  item_intercepts <- sprintf("y.%1$d ~ nu.%1$d*1;", t) |>
    paste(collapse = "\n")
  item_residuals <- sprintf("y.%1$d ~~ epsilon.%1$d*y.%1$d;", t) |>
    paste(collapse = "\n")

  # Define Latent Factors
  define_factors <- lapply(theta, function(x){
    paste0(x," =~ ",
           paste(sprintf(paste0(x,".%1$d*y.%1$d"), t),
                 collapse = " + "),";")
  }) |> paste(collapse = "\n")
  factor_means <- sprintf("%s ~ 0*1;", theta)|>
    paste(collapse = "\n")
  factor_covar <- outer(theta,theta, "paste", sep = " ~~ ")[
    lower.tri(outer(theta,theta, "paste", sep = " ~~ "), diag = TRUE)
  ] |> paste(collapse = "; \n")

  ## Define Phantoms for Latent Means
  define_phantoms <- lapply(theta, function(x){
    sprintf("%1$s_ph =~ 0; %1$s_ph ~ NA*1 + label(%1$s_m)*1 + %2$s?1;",
            x, start_values[x])
  }) |> paste(collapse = "\n")

  # Define Constraints
  ## Item Intercept Constraints
  nu_constraints <- sprintf(
    paste0("nu.%1$d == ", gsub("x_ti", "%2$f", mean_fx),";"), t, t-center_obs
  ) |> paste(collapse = "\n")

  ## Factor Loading Constraints
  lambda_constraints <- lapply(paste0(theta,"_m"), function(x){
    sprintf(paste0(
      substr(x, 1, nchar(x)-2),
      ".%1$d == ",
      D(mean_fx, x) |>
        deparse1(width.cutoff = 500, collapse = "") |>
        gsub("x_ti","%2$f",x = _),
      ";"
    ), t, t-center_obs) |> paste(collapse = "\n")
  }) |> paste(collapse = "\n\n")

  model.syntax <-paste(
    c("# Define Factors",define_factors, factor_means, factor_covar,
      "# Define Phantoms", define_phantoms,
      "# Define Items", item_intercepts, item_residuals,
      "# Define Constraints",nu_constraints, lambda_constraints),
    collapse = "\n\n")
  if(verbose) cat(model.syntax)
  return(model.syntax)
}

tidymessage <- function(..., prefix = " ", initial = ""){
  message(strwrap(..., prefix = prefix, initial = initial))
}
