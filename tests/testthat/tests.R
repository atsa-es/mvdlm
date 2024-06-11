
test_that("model fit with time varying formula works", {
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
                    data = data, chains = 1, iter = 20)

expect(class(fit$fit) == "stanfit")
})

test_that("model fit with normal, gamma, and lognormal family works", {
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
  time_varying = y ~ -1 + cov1 + cov2
  formula = y ~ 1

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 est_df = FALSE,
                 family = c("normal"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")

  data$y = data$y + abs(min(data$y)) + 1

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 est_df = FALSE,
                 family = c("gamma"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 est_df = FALSE,
                 family = c("lognormal"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")
})

test_that("model fit with binomial family works", {
  set.seed(123)
  N = 20
  data = data.frame("y" = runif(N),
                    "cov1" = rnorm(N),
                    "cov2" = rnorm(N),
                    "year" = 1:N,
                    "season" = sample(c("A","B"), size=N, replace=TRUE))
  b_1 = cumsum(rnorm(N))
  b_2 = cumsum(rnorm(N))
  data$y = rbinom(n=nrow(data), prob=plogis(data$cov1*b_1 + data$cov2*b_2), size=1)
  time_varying = y ~ -1 + cov1 + cov2
  formula = y ~ 1

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 est_df = FALSE,
                 family = c("binomial"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")
})

test_that("model fit with poisson and negative binomial family works", {
  set.seed(123)
  N = 20
  data = data.frame("y" = runif(N),
                    "cov1" = rnorm(N),
                    "cov2" = rnorm(N),
                    "year" = 1:N,
                    "season" = sample(c("A","B"), size=N, replace=TRUE))
  b_1 = cumsum(rnorm(N))
  b_2 = cumsum(rnorm(N))
  data$y = rpois(n=nrow(data), exp(data$cov1*b_1 + data$cov2*b_2))
  time_varying = y ~ -1 + cov1 + cov2
  formula = y ~ 1

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 est_df = FALSE,
                 family = c("nbinom2"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")
})

test_that("estimating df works", {
  set.seed(123)
  N = 20
  data = data.frame("y" = runif(N),
                    "cov1" = rnorm(N),
                    "cov2" = rnorm(N),
                    "year" = 1:N,
                    "season" = sample(c("A","B"), size=N, replace=TRUE))
  b_1 = cumsum(rnorm(N))
  b_2 = cumsum(rnorm(N))
  data$y = rpois(n=nrow(data), exp(data$cov1*b_1 + data$cov2*b_2))
  time_varying = y ~ -1 + cov1 + cov2
  formula = y ~ 1

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 est_df = TRUE,
                 family = c("nbinom2"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")
})


test_that("uncorrelated random walk works", {
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
  time_varying = y ~ -1 + cov1 + cov2
  formula = y ~ 1

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 correlated_rw = FALSE,
                 family = c("normal"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")
})

# test_that("plotting works", {
#   set.seed(123)
#   N = 20
#   data = data.frame("y" = runif(N),
#                     "cov1" = rnorm(N),
#                     "cov2" = rnorm(N),
#                     "year" = 1:N,
#                     "season" = sample(c("A","B"), size=N, replace=TRUE))
#   b_1 = cumsum(rnorm(N))
#   b_2 = cumsum(rnorm(N))
#   data$y = data$cov1*b_1 + data$cov2*b_2
#   time_varying = y ~ -1 + cov1 + cov2
#   formula = y ~ 1
#
#   fit <- fit_dlm(formula = formula,
#                  time_varying = time_varying,
#                  time = "year",
#                  family = c("normal"),
#                  data = data, chains = 1, iter = 20)
#
#   g <- dlm_trends(fit)
#   expect_equal(class(g), "list")
#   expect_equal(names(g)[1], "plot")
#   expect_equal(names(g)[2], "b_varying")
#   expect_equal(prod(dim(g[[2]])), 200)
# })

test_that("missing covariates works", {
  set.seed(123)
  N = 20
  data = data.frame("y" = runif(N),
                    "cov1" = rnorm(N),
                    "cov2" = rnorm(N),
                    "year" = 1:N,
                    "season" = sample(c("A","B"), size=N, replace=TRUE))
  b_1 = cumsum(rnorm(N))
  b_2 = cumsum(rnorm(N))
  data$y = rpois(n=nrow(data), exp(data$cov1*b_1 + data$cov2*b_2))
  data$cov1[10] = NA
  data$cov2[3] = NA
  time_varying = y ~ cov1 + cov2
  formula = NULL

  fit <- fit_dlm(formula = formula,
                 time_varying = time_varying,
                 time = "year",
                 family = c("normal"),
                 data = data, chains = 1, iter = 20)
  expect(class(fit$fit) == "stanfit")
})

test_that("model fit missing data row works", {
  set.seed(123)
  N = 20
  data = data.frame("y" = runif(N),
                    "cov1" = rnorm(N),
                    "cov2" = rnorm(N),
                    "year" = 1:N,
                    "season" = sample(c("A","B"), size=N, replace=TRUE))
  data = data[-18,]
  N = nrow(data)
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
                 data = data, chains = 1, iter = 20)

  expect(class(fit$fit) == "stanfit")
})


test_that("model fit missing observation works", {
  set.seed(123)
  N = 20
  data = data.frame("y" = runif(N),
                    "cov1" = rnorm(N),
                    "cov2" = rnorm(N),
                    "year" = 1:N,
                    "season" = sample(c("A","B"), size=N, replace=TRUE))
  data$y[18] = NA

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
                 data = data, chains = 1, iter = 20)

  expect(class(fit$fit) == "stanfit")
})

test_that("dlm_trends works", {
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
  #d <- try(dlm_trends(fit), silent=TRUE)
  d <- dlm_trends(fit)
  #expect(class(d) == "list")
  expect_equal(as.numeric(round(d$b_varying$estimate[1:5], 2)), c(0.31,0.12, -0.35, -0.69, -0.62))
})

