leni <- function(
    model = NULL,
    target_fx = "quadratic",
    theta = c("a0","ax","ay"),
    bootstrap = FALSE,
    model.class = "lme",
    data = NULL
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
    dat <- NULL
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

  # Identify Model Type - old
  fixed_effects <- leni_fixed_effects(model = model,
                                      model.class = model.class,
                                      ...)
  if(!is.null(modx)){
    cond_effects <- leni_conditional_effects(model = model,
                                             model.class = model.class,
                                             moderator = modx,
                                             ...)
  }
  if(model.type == "lme"){
    leni_random_effects(model = model, ...)
  }
}
