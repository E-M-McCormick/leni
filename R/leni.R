target_fx = "cubic"; theta = c("xN","yN","d","h"); bootstrap = TRUE; model.class = "lm"; data = NULL
leni <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = FALSE,
    model.class = "lme",
    data = NULL,
    coef.idx = NULL,
    modx.idx = NULL,
    ...
){

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
  }
  if(is.null(data)){
    bootstrap <- FALSE
    tidymessage("The raw data is required to obtain bootstrap results.
                 Defaulting to analytic standard errors.")
  }

  # Try to automatically parse target function
  read.coefs <- if(model.class == "lm"){
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

  # Obtain Model Estimates
  fixed_effects <- leni_fixed_effects(model = model,
                                      model.class = model.class,
                                      ...)
  if(!is.null(modx)){
    cond_effects <- leni_conditional_effects(model = model,
                                             model.class = model.class,
                                             moderator = modx,
                                             ...)
  }
  if(model.class == "lme"){
    leni_random_effects(model = model, ...)
  }
}
