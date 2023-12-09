leni <- function(
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
    varcov = TRUE,
    ci = 0.95,
    verbose = FALSE,
    ...
){

  # Generate Output Structure
  output <- list()
  output[['function_call']] <- as.list(sys.frame(which = 1))
  output[['function_call']][['model']] <- NULL
  output[['linear_model']] <- model

  # Default Arguments
  if(is.null(model)){stop(tidymessage("No model object provided..."))}
  if(is.logical(bootstrap) & bootstrap){bootstrap <- 5000}

  # Identify Model Type
  if(is.null(model.class)){model.class <- class(model)}
  if(grepl(model.class, "lm")){
    model.class <- "lm"
  } else if(any(grepl(model.class[1],
                      c("lmerMod","lmerModLmerTest","nlme","lme")))){
    model.class <- "lme"
  }

  # Extract Raw Data for Bootstrap
  if(is.numeric(bootstrap) & bootstrap > 0){
    if(grepl(model.class, "lm")){data <- model$model}
    if(class(model) %in% c("lmerMod","lmerModLmerTest")){data <- model@frame}
    if(class(model) %in% c("nlme","lme")){data <- nlme::getData(model)}
    if(is.null(data)){
      bootstrap <- FALSE
      tidymessage("The raw data is required to obtain bootstrap results.
                 Defaulting to analytic standard errors.")
    }
  }

  # Try to automatically parse target function
  read.coefs <<- if(model.class == "lm"){
    stats::coef
  } else if(model.class == "lme"){
    lme4::fixef
  }

  if(is.null(target_fx)){
    if(length(read.coefs(model)) == 3) target_fx <- "quadratic"
    if(any(grepl("\\^3)",names(read.coefs(model)))) |
       any(grepl("poly(, 3)*",names(read.coefs(model))))) {
      target_fx <- "cubic"
    } else if(any(grepl("\\^2)",names(read.coefs(model)))) |
              any(grepl("poly(, 2)*",names(read.coefs(model))))){
      target_fx <- "quadratic"
    }
  }
  if(is.null(target_fx)){stop(tidymessage('Please specify target function
                                          (i.e., "quadratic" or "cubic").'))}

  # Specify Coefficient Index
  if (is.null(coef.idx)){
    if (grepl(target_fx, "quadratic")){
      if (verbose) tidymessage("Defaulting to model coefficients 1-3
                               as fixed effects targets. ")
      coef.idx <- c(1:3)
    } else if (grepl(target_fx, "cubic")){
      if (verbose) tidymessage("Defaulting to model coefficients 1-4
                               as fixed effects targets.")
      coef.idx <- c(1:4)
    }
  }

  # Specify Moderator Index
  if (!is.null(modx) & is.null(modx.idx)){
    if (grepl(target_fx, "quadratic")){
      if (verbose) tidymessage('Defaulting to first 3 coefficients containing
                                the moderator as targets. Specify "mod.idx" if
                                an alternative is needed.')
      modx.idx <- grep(modx, names(read.coefs(model)))[1:3]
    } else if (grepl(target_fx, "cubic")){
      if (verbose) tidymessage('Defaulting to first 4 coefficients containing
                                the moderator as targets. Specify "mod.idx" if
                                an alternative is needed.')
      modx.idx <- grep(modx, names(read.coefs(model)))[1:4]
    }
  }

  # Obtain Model Estimates
  if(is.null(modx)){
    output[['leni_model']][['fixed_effects']] <-
      leni_fixed_effects(model = model,
                         target_fx = target_fx,
                         theta = theta,
                         bootstrap = bootstrap,
                         bootSeed = bootSeed,
                         model.class = model.class,
                         data = data,
                         coef.idx = coef.idx,
                         ci = ci,
                         ...)
  }
  if(!is.null(modx)){
    output[['leni_model']][['fixed_effects']] <-
      leni_conditional_effects(model = model,
                               target_fx = target_fx,
                               theta = theta,
                               bootstrap = bootstrap,
                               bootSeed = bootSeed,
                               model.class = model.class,
                               data = data,
                               coef.idx = coef.idx,
                               modx = modx,
                               modx.idx = modx.idx,
                               ci = ci,
                               ...)
  }
  random_effects <- NULL
  if(model.class == "lme" & varcov){
    output[['leni_model']][['random_effects']] <-
      leni_random_effects(model = model,
                          target_fx = target_fx,
                          theta = theta,
                          bootstrap = bootstrap,
                          bootSeed = bootSeed,
                          model.class = model.class,
                          data = data,
                          coef.idx = coef.idx,
                          ci = ci,
                          ...)
  }

  return(output)
}
