leni_sem <- function(
    target_fx = "cubic",
    theta = NULL,
    start_values = NULL,
    number_obs = 5,
    center_obs = 1,
    spacing = 1,
    verbose = FALSE,
    y.name = NULL,
    time.name = NULL
){
  # Check Initial arguments
  if(is.null(target_fx)){stop('Target function is not defined. Stopping...')}
  if(is.null(y.name)){y.name = "y."}
  if(is.null(time.name)){time.name = "x_ti"}
  if(any(grepl(time.name, theta))){stop('`time.name` argument is a substring
                                         of at least one element in `theta`,
                                         which will cause weird results.
                                         Please alter before proceeding.')}

  # Define Time Points
  t <- if(length(spacing) == 1){
    seq(1, (spacing*number_obs), by = spacing)
  } else if(length(spacing) > 1) {
    c(1, 1 + cumsum(spacing))
  }

  # Define Target Function & Parameters
  c(fx, theta) %tin% nonlinear_function_library(target_fx = target_fx,
                                                theta = theta,
                                                time.name = time.name)
  if(!grepl(time.name, fx)){
    stop(tidymessage('Unable to find the time-related variable within the target
                      function (e.g., "x_ti"). Please use the time.name
                      argument to specify if you are using a custom
                      variable name.'))}

  # Define Start Values
  if(is.null(start_values)){start_values <- rep(0, length(theta))
  } else if(length(start_values) != length(theta)){
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
  item_intercepts <- sprintf("%1$s%2$g ~ nu.%2$g*1;", y.name, t) |>
    paste(collapse = "\n")
  item_residuals <- sprintf("%1$s%2$g ~~ epsilon.%2$g*%1$s%2$g;", y.name, t) |>
    paste(collapse = "\n")

  # Define Latent Factors
  define_factors <- lapply(theta, function(x){
    paste0(x," =~ ",
           paste(sprintf(paste0(x,".%2$g*%1$s%2$g"), y.name, t),
                 collapse = " + "),";")
  }) |> paste(collapse = "\n")
  factor_means <- sprintf("%s ~ 0*1;", theta)|>
    paste(collapse = "\n")
  factor_covar <- paste0(outer(theta,theta, "paste", sep = " ~~ ")[
    lower.tri(outer(theta,theta, "paste", sep = " ~~ "), diag = TRUE)
  ] |> paste(collapse = "; \n"),";")

  ## Define Phantoms for Latent Means
  define_phantoms <- lapply(theta, function(x){
    sprintf('%1$s_ph =~ 0; %1$s_ph ~ NA*1 + label("%1$s_m")*1 + %2$s?1;',
            x, start_values[x])
  }) |> paste(collapse = "\n")

  # Define Constraints
  ## Item Intercept Constraints
  nu_constraints <- sprintf(
    paste0("nu.%1$g == ", gsub(time.name, "%2$g", mean_fx),";"), t, t-center_obs
  ) |> paste(collapse = "\n")

  ## Factor Loading Constraints
  lambda_constraints <- lapply(paste0(theta,"_m"), function(x){
    if(grepl(time.name,
             stats::D(mean_fx, x) |>
             deparse1(width.cutoff = 500, collapse = ""),
             fixed = TRUE)){
      sprintf(paste0(
        substr(x, 1, nchar(x)-2),
        ".%1$g == ",
        stats::D(mean_fx, x) |>
          deparse1(width.cutoff = 500, collapse = "") |>
          gsub(time.name,"%2$g",x = _),
        ";"
      ), t, t-center_obs) |> paste(collapse = "\n")
    } else {
      sprintf(paste0(
        substr(x, 1, nchar(x)-2),
        ".%1$g == ",
        stats::D(mean_fx, x) |>
          deparse1(width.cutoff = 500, collapse = ""),
        ";"
      ), t) |> paste(collapse = "\n")
    }
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

