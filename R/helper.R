# Helper Functions
`%tin%` <- function(x, y) {
  mapply(assign, as.character(substitute(x)[-1]), y,
         MoreArgs = list(envir = parent.frame()))
  invisible()
}

tidymessage <- function(..., prefix = " ", initial = ""){
  message(strwrap(..., prefix = prefix, initial = initial))
}

citations <- list(
  CdT_2002 = 'Cudeck, R., & Du Toit, S. H. C. (2002).
  A version of quadratic regression with interpretable parameters.
  Multivariate Behavioral Research, 37(4), 501–519.',
  Mcc_2023 = 'McCormick (2023).
  Deriving models of change with interpretable parameters:
  linear estimation with nonlinear inference. PsyArixv.
  htts://doi.org'
)

nonlinear_function_library <- function(target_fx = NULL, theta = NULL,
                                       time.name = NULL){
  ## Quadratic Functions
  if(grepl(target_fx,"quadratic")){
    tidymessage('Defaulting to the standard alternative quadratic from
                 Cudeck & du Toit (2002). For the other alternatives
                 please use: target_fx = "quadratic_gamma_a0",
                "quadratic_gamma_ay", or "quadratic_ac"')
    fx <- paste0("ay - (ay - a0)*((",time.name," / ax) - 1)^2") |>
      str2expression()
    theta <- c("a0","ax","ay")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["CdT_2002"]])

  } else if(grepl(target_fx,"quadratic_gamma_a0")) {
    fx <- paste0("a0 - gamma*(((",time.name," / ax) - 1)^2 - 1)") |>
      str2expression()
    theta <- c("a0","ax","gamma")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["CdT_2002"]])

  } else if(grepl(target_fx,"quadratic_gamma_ay")){
    fx <- paste0("ay - gamma*((",time.name," / ax) - 1)^2") |>
      str2expression()
    theta <- c("ax","ay","gamma")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["CdT_2002"]])

  } else if(grepl(target_fx,"quadratic_ac")){
    fx <- paste0("ay + ac*(",time.name," - ax)^2") |>
      str2expression()
    theta <- c("ax","ay","ac")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["Mcc_2023"]])


  ## Cubic Functions
  } else if(grepl(target_fx,"cubic") |
            grepl(target_fx, "cubic_h")){
    tidymessage('Defaulting to the alternative cubic with h from
                 McCormick (2024). For the version with betaN
                 please use: target_fx = "cubic_betaN"')
    fx <- paste0("yN - (h/2)*(((",time.name,"-xN)/d)^3 -
                 3*((",time.name,"-xN)/d))") |>
      str2expression()
    theta <- c("xN","yN","d","h")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["Mcc_2023"]])

  } else if(grepl(target_fx,"cubic_betaN")){
    fx <- paste0("yN - (betaN*d/3)*(((",time.name,"-xN)/d)^3 -
                     3*((",time.name,"-xN)/d))") |>
      str2expression()
    theta <- c("xN","yN","d","betaN")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["Mcc_2023"]])

  ## Multiphase Functions
  } else if(grepl(target_fx, "cubic_multiphase") |
            grepl(target_fx, "cubic_multiphase_h")){
    tidymessage('Defaulting to the alternative multiphase cubic with h from
                 McCormick (2024). For the version with betaN
                 please use: target_fx = "cubic_multiphase_betaN"')
    fx <- paste0(
      "yN - (h/2)*(
        (((1/2)*(sqrt((xN-d-",time.name,")^2)-
        sqrt((",time.name,"-xN-d)^2)))/d)^3 -
          3*(((1/2)*(sqrt((xN-d-",time.name,")^2)-
      sqrt((",time.name,"-xN-d)^2)))/d))"
    ) |> str2expression()
    theta <- c("xN","yN","d","h")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["Mcc_2023"]])

  } else if(grepl(target_fx, "cubic_multiphase_betaN")){
    fx <- paste0(
      "yN - (betaN*d/3)*(
        (((1/2)*(sqrt((xN-d-",time.name,")^2)-
        sqrt((",time.name,"-xN-d)^2)))/d)^3 -
          3*(((1/2)*(sqrt((xN-d-",time.name,")^2)-
      sqrt((",time.name,"-xN-d)^2)))/d))"
    ) |> str2expression()
    theta <- c("xN","yN","d","betaN")
    tidymessage('Please considering citing the following:')
    tidymessage(citations[["Mcc_2023"]])

  ## Logistic Functions
  } else if(grepl(target_fx, "1PL")){
    fx <- paste0("(1/(1+(",time.name,"/xN)^(-xN/2)))") |>
      str2expression()
    theta <- c("xN")
  } else if(grepl(target_fx, "2PL")){
    fx <- paste0("(1/(1+(",time.name,"/xN)^(-hill)))") |>
      str2expression()
    theta <- c("xN","hill")
  } else if(grepl(target_fx, "3PL")){
    fx <- paste0("lower+((1-lower)/(1+(",time.name,"/xN)^(-hill)))") |>
      str2expression()
    theta <- c("lower","xN","hill")
  } else if(grepl(target_fx, "4PL")){
    fx <- paste0("lower+((upper-lower)/(1+(",time.name,"/xN)^(-hill)))") |>
      str2expression()
    theta <- c("lower","upper","xN","hill")


  ## Stop on Improper target_fx inputs
  } else if(!is.character(target_fx) & !is.expression(target_fx)) {
    stop(tidymessage('Unable to parse target function as a
                     string for conversion to expression.'))

  ## Custom Target Function
  } else {
    if(is.null(theta)){stop(tidymessage(
      "Parameters not specified in theta. Unable to parse target function."
    ))}
    fx <- str2expression(target_fx)
    theta <- theta
  }
  return(list(fx, theta))
}
