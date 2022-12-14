% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlm_trends.R
\name{dlm_trends}
\alias{dlm_trends}
\title{Summarize and plot time varying coefficients from the fitted model}
\usage{
dlm_trends(fitted_model)
}
\arguments{
\item{fitted_model}{A fitted model object}
}
\value{
A list containing the plot and data used
to fit the model. These include \code{plot} and \code{b_varying}
}
\description{
Summarize and plot time varying coefficients from the fitted model
}
\examples{
\donttest{
set.seed(123)
N = 20
data = data.frame("y" = runif(N),
                  "cov1" = rnorm(N),
                  "cov2" = rnorm(N),
                  "year" = 1:N,
                  "season" = sample(c("A","B"), size=N, replace=TRUE))
b_1 = cumsum(rnorm(N))
b_2 = cumsum(rnorm(N))
data$y = data$cov1*b_1 + data$cov2*b_2
time_varying = y ~ cov1 + cov2
formula = NULL
fit <- fit_dlm(formula = formula,
               time_varying = time_varying,
               time = "year",
               est_df = FALSE,
               family = c("normal"),
               data=data, chains = 1, iter = 20)
dlm_trends(fit)
}

}
