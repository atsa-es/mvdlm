#' Fit a Bayesian multivariate dynamic linear model with Stan
#'
#' Fit a Bayesian multivariate dynamic linear model with Stan that optionally includes covariates to estimate
#' effects, extremes (Student-t distribution), etc.
#'
#' @param formula The model formula for the fixed effects; at least this formula or `time_varying` needs to have the response included
#' @param time_varying The model formula for the time-varying effects; at least this formula or `formula` needs to have the response included
#' @param time String describing the name of the variable corresponding to time, defaults to "year"
#' @param est_df Whether or not to estimate deviations of B as Student - t with estimated degrees of freedom, defaults to `FALSE`
#' @param family, The name of the family used for the response; can be one of "normal","binomial","possion","nbinom2","gamma","lognormal"
#' @param correlated_rw, Whether to estimate time-varying parameters as correlated random walk, defaults to TRUE
#' @param data The data frame including response and covariates for all model components
#' @param chains Number of mcmc chains, defaults to 3
#' @param iter Number of mcmc iterations, defaults to 2000
#' @param warmup Number iterations for mcmc warmup, defaults to 1/2 of the iterations
#' @param ... Any other arguments to pass to [rstan::sampling()].
#' @export
#' @return A list containing the fitted model and arguments and data used
#' to fit the model. These include `model` (the fitted model object of class `stanfit`),

#' @importFrom rstan sampling
#' @importFrom stats model.frame model.matrix model.response
#' @import Rcpp
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' N = 20
#' data = data.frame("y" = runif(N),
#'                   "cov1" = rnorm(N),
#'                   "cov2" = rnorm(N),
#'                   "year" = 1:N,
#'                   "season" = sample(c("A","B"), size=N, replace=TRUE))
#' b_1 = cumsum(rnorm(N))
#' b_2 = cumsum(rnorm(N))
#' data$y = data$cov1*b_1 + data$cov2*b_2
#' time_varying = y ~ cov1 + cov2
#' formula = NULL
#'
#' # fit a model with a time varying component
#' fit <- fit_dlm(formula = formula,
#'                time_varying = time_varying,
#'                time = "year",
#'                est_df = FALSE,
#'                family = c("normal"),
#'                data=data, chains = 1, iter = 20)
#'
#' # fit a model with a time varying and fixed component (here, fixed intercept)
#' fit <- fit_dlm(formula = y ~ 1,
#'                time_varying = y ~ -1 + cov1 + cov2,
#'                time = "year",
#'                est_df = FALSE,
#'                family = c("normal"),
#'                data=data, chains = 1, iter = 20)
#'
#' #' # fit a model with deviations modeled with a multivariate Student-t
#' fit <- fit_dlm(formula = y ~ 1,
#'                time_varying = y ~ -1 + cov1 + cov2,
#'                time = "year",
#'                est_df = TRUE,
#'                family = c("normal"),
#'                data=data, chains = 1, iter = 20)
#'
#' #' #' # fit a model with deviations modeled with a multivariate Student-t
#' fit <- fit_dlm(formula = y ~ 1,
#'                time_varying = y ~ -1 + cov1 + cov2,
#'                time = "year",
#'                est_df = TRUE,
#'                family = c("normal"),
#'                data=data, chains = 1, iter = 20)
#' }
#'
fit_dlm <- function(formula = NULL,
                    time_varying = NULL,
                    time = "year",
                    est_df = FALSE,
                    family = c("normal", "binomial", "poisson", "nbinom2", "gamma", "lognormal"),
                    correlated_rw = TRUE,
                    data,
                    chains = 3,
                    iter = 2000,
                    warmup = floor(iter / 2),
                    ...) {


  # first fill in any missing years
  all_years <- seq(min(data[[time]]), max(data[[time]]))
  missing_yrs <- all_years[which(all_years %in% data[[time]] == FALSE)]
  if(length(missing_yrs) > 0) {
    df_new = data[1:length(missing_yrs),]
    for(i in 1:ncol(df_new)) df_new[,i] = NA
    df_new[[time]] = missing_yrs
    data = rbind(data, df_new)
  }
  # add intercept column to data
  data$`(Intercept)` <- 1

  recognized_families <- c("normal", "binomial", "poisson", "nbinom2", "gamma", "lognormal")
  family <- family[1]
  if (family %in% recognized_families == FALSE) {
    stop("Error: family not recognized")
  } else {
    family <- match(family, recognized_families)
  }

  # parse formulas
  est_fixed_coef <- FALSE
  est_varying_coef <- FALSE
  n_fixed <- 0
  n_varying <- 0

  y <- NULL
  tv_pars <- NULL
  fixed_pars <- NULL
  if (!is.null(formula)) {
    model_frame <- model.frame(formula, data, na.action=na.pass)
    y <- model.response(model_frame)
    model_matrix <- model.matrix(formula, model_frame)
    fixed_pars <- colnames(model_matrix)
    est_fixed_coef <- TRUE
    fixed_dat <- cbind(model_matrix, c(data[, time]))
    colnames(fixed_dat)[ncol(fixed_dat)] <- "time"
    fixed_dat[,ncol(fixed_dat)] = fixed_dat[,ncol(fixed_dat)] - min(fixed_dat[,ncol(fixed_dat)]) + 1
    n_fixed <- ncol(fixed_dat) - 1
    fixed_time <- rep(fixed_dat[, "time"], ncol(fixed_dat) - 1)
    fixed_var <- sort(rep(1:n_fixed, nrow(fixed_dat)))
    fixed_x <- c(as.matrix(fixed_dat[, which(colnames(fixed_dat) != "time")]))
    fixed_N <- length(fixed_time)
    n_fixed_NAs <- length(which(is.na(fixed_x)))
    fixed_NAs <- 0 # dummy
    if (n_fixed_NAs > 0) {
      fixed_NAs <- c(which(is.na(fixed_x)), 0, 0)
      fixed_x[which(is.na(fixed_x))] = 0
    } else {
      fixed_NAs <- c(0, 0)
    }
  } else {
    n_fixed <- 0
    fixed_time <- c(0, 0)
    fixed_var <- c(0, 0)
    fixed_x <- c(0, 0)
    fixed_N <- 2
    n_fixed_NAs <- 0
    fixed_NAs <- c(0, 0) # dummy
  }
  if (!is.null(time_varying)) {
    model_frame <- model.frame(time_varying, data, na.action=na.pass)
    if (is.null(y)) y <- model.response(model_frame)
    model_matrix <- model.matrix(time_varying, model_frame)
    tv_pars <- colnames(model_matrix)
    est_varying_coef <- TRUE
    varying_dat <- cbind(model_matrix, c(data[, time]))
    colnames(varying_dat)[ncol(varying_dat)] <- "time"
    varying_dat[,ncol(varying_dat)] <- varying_dat[,ncol(varying_dat)] - min(varying_dat[,ncol(varying_dat)]) + 1
    n_varying <- ncol(varying_dat) - 1
    varying_time <- rep(varying_dat[, "time"], ncol(varying_dat) - 1)
    varying_var <- sort(rep(1:n_varying, nrow(varying_dat)))
    varying_x <- c(as.matrix(varying_dat[, which(colnames(varying_dat) != "time")]))
    varying_N <- length(varying_time)
    n_varying_NAs <- length(which(is.na(varying_x)))
    varying_NAs <- 0 # dummy
    if (n_varying_NAs > 0) {
      varying_NAs <- c(which(is.na(varying_x)), 0, 0)
      varying_x[which(is.na(varying_x))] = 0
    } else {
      varying_NAs <- c(0, 0)
    }
  } else {
    n_varying <- 0
    varying_time <- 0
    varying_var <- 0
    varying_x <- 0
    varying_N <- 1
    n_varying_NAs <- 0
    varying_NAs <- c(0, 0) # dummy
  }

  stan_data <- list(
    y = y[which(!is.na(y))],
    y_int = as.integer(y[which(!is.na(y))]),
    N = length(which(!is.na(y))),
    y_indx = which(!is.na(y)),
    nT = max(c(fixed_time, varying_time), na.rm = T),
    est_fixed = as.numeric(est_fixed_coef),
    est_varying = as.numeric(est_varying_coef),
    n_fixed_covars = n_fixed,
    fixed_N = fixed_N,
    fixed_time_indx = fixed_time,
    fixed_var_indx = fixed_var,
    fixed_x_value = fixed_x,
    n_varying_covars = n_varying,
    varying_N = varying_N,
    varying_time_indx = varying_time,
    varying_var_indx = varying_var,
    varying_x_value = varying_x,
    est_df = as.numeric(est_df),
    family = family,
    n_fixed_NAs = n_fixed_NAs,
    fixed_NAs = fixed_NAs,
    n_varying_NAs = n_varying_NAs,
    varying_NAs = varying_NAs,
    correlated_rw = as.numeric(correlated_rw)
  )

  pars <- c("eta", "sigma", "log_lik", "lp__")
  if(est_varying_coef == TRUE) pars <- c(pars, "b_varying")
  if(est_fixed_coef == TRUE) pars <- c(pars, "b_fixed")
  if(family %in% c("normal","negbin2","gamma","lognormal")) pars <- c(pars, "phi")
  if(est_df == TRUE) pars <- c(pars, "nu")
  if(correlated_rw == TRUE) pars <- c(pars, "R", "Sigma", "Lcorr")

  sampling_args <- list(
    object = stanmodels$dlm,
    chains = chains,
    iter = iter,
    warmup = warmup,
    pars = pars,
    data = stan_data, ...
  )
  fit <- do.call(sampling, sampling_args)

  return(list(
    fit = fit,
    "fixed_pars" = fixed_pars,
    "time_varying_pars" = tv_pars,
    fixed_formula = formula,
    time_varying_formula = time_varying,
    time = time,
    est_df = est_df,
    stan_data = stan_data,
    raw_data = data
  ))
}
