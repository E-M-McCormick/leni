organize_output <- function(
    model,
    model.class = NULL,
    temp_call = NULL,
    fixed_effects = NULL,
    random_effects = NULL,
    bootstrap = FALSE,
    modx = NULL,
    ci = 0.95
){
  if(is.null(fixed_effects)) stop(tidymessage("No fixed effects output detected.
                                               Stopping..."))


  if (grepl(model.class, "lm")){

  } else if (grepl(model.class, "lme")){
    crit.val <- qt((1 - ci) / 2,
                   df = length(unique(leni_realdata@flist[[1]])),
                   lower.tail = FALSE)

    if (is.numeric(bootstrap) & bootstrap > 0){
      bootstrap_estimates <- data.frame(
        if(!is.null(modx)){
          sub("modx", modx, names(fixed_effects$fixed_effects))
        } else {
          names(fixed_effects$fixed_effects)
        },
        est = fixed_effects$fixed_effects,
        se = fixed_effects$fixed_effects_se,
        ci.lower = sapply(1:ncol(fixed_effects$bootstrap_samples$t),
                          function(x){
                            boot::boot.ci(temp$fixed_effects$bootstrap_samples,
                                          type = "bca",
                                          index = x)$bca[4]
                          }),
        ci.upper = sapply(1:ncol(fixed_effects$bootstrap_samples$t),
                          function(x){
                            boot::boot.ci(temp$fixed_effects$bootstrap_samples,
                                          type = "bca",
                                          index = x)$bca[5]
                          }),
        est.r = fixed_effects$fixed_effects_robust,
        se.r = fixed_effects$fixed_effects_se_robust
      )
    } else {
      estimates <- data.frame(
        row.names = if(!is.null(modx)){
          sub("modx", modx, names(fixed_effects$fixed_effects))
        } else {
          names(fixed_effects$fixed_effects)
        },
        est = fixed_effects$fixed_effects,
        se = fixed_effects$fixed_effects_se,
        ci.lower = fixed_effects$fixed_effects -
          crit.val*fixed_effects$fixed_effects_se,
        ci.upper = fixed_effects$fixed_effects +
          crit.val*fixed_effects$fixed_effects_se
      )
    }


  }


}
