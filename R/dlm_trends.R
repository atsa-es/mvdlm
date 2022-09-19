#' Summarize and plot time varying coefficients from the fitted model
#'
#' @param fitted_model A fitted model object
#' @export
#' @return A list containing the plot and data used
#' to fit the model. These include `plot` and `b_varying`

#' @importFrom broom.mixed tidy
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' N = 20
#' data = data.frame("y" = runif(N),
#'                   "cov1" = rnorm(N),
#'                   "cov2" = rnorm(N),
#'                   "year" = 1:N,
#'                   "season" = sample(c("A","B"), size=N, replace=T))
#' b_1 = cumsum(rnorm(N))
#' b_2 = cumsum(rnorm(N))
#' data$y = data$cov1*b_1 + data$cov2*b_2
#' time_varying = y ~ cov1 + cov2
#' formula = NULL
#' fit <- fit_dlm(formula = formula,
#'                time_varying = time_varying,
#'                time = "year",
#'                est_df = FALSE,
#'                family = c("normal"),
#'                data, chains = 1, iter = 20)
#' dlm_trends(fit)
#' }
#'
dlm_trends <- function(fitted_model) {

  tidy_pars <- broom.mixed::tidy(fitted_model$fit)

  indx <- grep("b_varying", tidy_pars$term)
  if(length(indx) == 0) {
    stop("Error: time varying parameters not found")
  }

  b_varying = tidy_pars[indx,] # subset
  b_varying$par <- rep(fit$time_varying_pars, each = fit$stan_data$nT) # add names
  b_varying$time <- rep(1:fit$stan_data$nT, length(fit$time_varying_pars))

  cols <- "#440154FF" # viridis::viridis(1)
  g <- ggplot(b_varying, aes(time, estimate)) +
    geom_ribbon(aes(ymin=estimate-1.96*std.error, ymax=estimate+1.96*std.error), fill=cols, alpha=0.5) +
    geom_line(col = cols) +
    facet_wrap(~par, scales="free_y") +
    ylab("Estimate") +
    xlab("Time") +
    theme_bw() +
    theme(strip.background =element_rect(fill="white"))
  return(list(plot = g, b_varying = b_varying))
}
