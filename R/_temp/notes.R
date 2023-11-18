####### Fixed Effects + Standard Errors
# Quadratic Model
linear.fit <- lm(y ~ 1 + I(x) + I(x^2), data = quadratic_ols)
ACOV_beta <- vcov(linear.fit)

f_a0 <- expression(b0)
f_ax <- expression(-b1/(2 * b2))
f_ay <- expression(b0 - (b1^2/(4 * b2)))
f_g <- expression(- (b1^2/(4 * b2)))
f_ac <- expression(b2)
c(b0,b1,b2) %tin% coef(linear.fit)

J <- matrix(  # see if I can simplify this J matrix code
  c(sapply(c("b0","b1","b2"), function(x) eval(D(f_a0,x))),
    sapply(c("b0","b1","b2"), function(x) eval(D(f_ax,x))),
    sapply(c("b0","b1","b2"), function(x) eval(D(f_ay,x))),
    sapply(c("b0","b1","b2"), function(x) eval(D(f_g,x))),
    sapply(c("b0","b1","b2"), function(x) eval(D(f_ac,x)))),
  nrow = 3, ncol = 5, byrow = FALSE) # check if dimensions are correct

## If Analytic Standard Errors
ACOV_alpha <- t(J) %*% ACOV_beta %*% J      # Find a way to slice J based on desired theta vector

## If Bootstrap Standard Errors

quad_boot <- function(data, sample, formula){
  fit <- lm(formula, data=data[sample,]) # change model fit for lm versus lme
  c(b0,b1,b2) %tin% coef(fit)[index] # slice index from the model object or user options

  return(setNames(
    c(b0, -b1/(2*b2), b0 - (b1^2)/(4*b2)),
    c("alpha_0", "alpha_x", "alpha_y"))) # slice this depending on desired theta
}

results <- boot::boot(data = dat, # from model object
                      statistic = quad_boot,
                      R = bootstrap,
                      formula = as.formula("y ~ 1 + x + x2")) # insert model equation from the model object

# For bootstrap option, return results$t0 and full results objects



-----------------------
# Cubic Model
f_xn <- expression(-b2 / (3 * b3))
f_yn <- expression(b0 - (b1 * b2) / (3 * b3) + (2 * b2^3) / (27 * b3^2))
f_d <- expression(sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2)))
f_h <- expression(-2 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^3)
f_bn <- expression(-3 * b3 * sqrt((b2^2 - 3 * b3 * b1) / (9 * b3^2))^2)
c(b0,b1,b2,b3) %tin% coef(linear.fit)

J <- matrix(
  c(sapply(c("b0","b1","b2","b3"), function(x) eval(D(f_xn, x))),
    sapply(c("b0","b1","b2","b3"), function(x) eval(D(f_yn, x))),
    sapply(c("b0","b1","b2","b3"), function(x) eval(D(f_d, x))),
    sapply(c("b0","b1","b2","b3"), function(x) eval(D(f_h, x))),
    sapply(c("b0","b1","b2","b3"), function(x) eval(D(f_bn, x)))),
  nrow = 4, ncol = 5, byrow = FALSE)

ACOV_theta <- t(J) %*% ACOV_beta %*% J      ## Find a way to slice J based on desired theta vector
