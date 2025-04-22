data {
  int N; // number of samples
  int y_indx[N]; // observations
  real y[N]; // observations
  int y_int[N]; // observations
  int nT; // number of time steps
  int fixed_N;
  int n_fixed_covars; // number fixed covariates
  int fixed_time_indx[fixed_N]; // time index of fixed effects
  int fixed_var_indx[fixed_N]; // variable index of fixed effects
  real fixed_x_value[fixed_N]; // value of covariates
  int varying_N;
  int n_varying_covars; // number fixed covariates
  int varying_time_indx[varying_N]; // time index of fixed effects
  int varying_var_indx[varying_N]; // variable index of fixed effects
  real varying_x_value[varying_N]; // value of covariates
  int est_df; // whether to estimate Student - t df
  int family;
  int n_fixed_NAs;
  int fixed_NAs[n_fixed_NAs + 2];
  int n_varying_NAs;
  int varying_NAs[n_varying_NAs + 2];
  int correlated_rw;
  real phi_scale; // scale for phi prior (e.g. SD for normal)
}
transformed data {
  vector[n_varying_covars] zeros;
  for(i in 1:n_varying_covars) zeros[i] = 0;
}
parameters {
  vector[n_fixed_covars] b_fixed;
  vector[n_varying_covars] b_devs0;
  vector[n_varying_covars] z_devs[nT - 1];  // ← added
  cholesky_factor_corr[n_varying_covars * correlated_rw] Lcorr;
  vector<lower=0>[n_varying_covars] sigma;
  vector<lower=0>[1] phi;
  vector<lower=0>[est_df] nu;
  vector[n_fixed_NAs] missing_fixed;
  vector[n_varying_NAs] missing_varying;
}
transformed parameters {
  vector[n_varying_covars] b_devs[nT - 1];
  vector[n_varying_covars] b_varying[nT];
  matrix[n_varying_covars * correlated_rw, n_varying_covars * correlated_rw] R;
  matrix[n_varying_covars * correlated_rw, n_varying_covars * correlated_rw] Sigma;
  vector[nT] eta;
  matrix[nT, n_fixed_covars] X_fixed;
  matrix[nT, n_varying_covars] X_varying;

  // Correlation structure
  if (correlated_rw == 1) {
    R = multiply_lower_tri_self_transpose(Lcorr);
    Sigma = quad_form_diag(R, sigma);
  }

  // Reassemble design matrices
  if (n_fixed_covars > 0) {
    for (i in 1:fixed_N) {
      X_fixed[fixed_time_indx[i], fixed_var_indx[i]] = fixed_x_value[i];
    }
    for (i in 1:n_fixed_NAs) {
      X_fixed[fixed_time_indx[fixed_NAs[i]], fixed_var_indx[fixed_NAs[i]]] = missing_fixed[i];
    }
  }

  if (n_varying_covars > 0) {
    for (i in 1:varying_N) {
      X_varying[varying_time_indx[i], varying_var_indx[i]] = varying_x_value[i];
    }
    for (i in 1:n_varying_NAs) {
      X_varying[varying_time_indx[varying_NAs[i]], varying_var_indx[varying_NAs[i]]] = missing_varying[i];
    }
  }

  // Non-centered parameterization for b_devs
  for (t in 1:(nT - 1)) {
    if (correlated_rw == 1) {
      // Correlated innovations
      b_devs[t] = Lcorr * (sigma .* z_devs[t]);
    } else {
      // Independent innovations
      b_devs[t] = sigma .* z_devs[t];
    }
  }

  // Time-varying coefficients from innovations
  b_varying[1] = b_devs0;
  for (t in 2:nT) {
    b_varying[t] = b_varying[t - 1] + b_devs[t - 1];
  }

  // Linear predictor
  for (t in 1:nT) eta[t] = 0;
  if (n_fixed_covars > 0) eta += X_fixed * b_fixed;
  if (n_varying_covars > 0) {
    for (t in 1:nT) {
      eta[t] += dot_product(X_varying[t], b_varying[t]);
    }
  }
}
model {
  // Priors
  sigma ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2.0);
  phi ~ normal(0, phi_scale);  // tighter obs noise prior
  b_fixed ~ normal(0, 1);
  nu ~ student_t(3, 0, 1);

  missing_fixed ~ normal(0, 1);
  missing_varying ~ normal(0, 1);
  b_devs0 ~ normal(0, 0.1);  // informative prior to anchor trajectory

  // Non-centered process priors
  for (t in 1:(nT - 1)) {
    if (est_df == 0) {
      to_vector(z_devs[t]) ~ normal(0, 1);
    } else {
      for (k in 1:n_varying_covars) {
        z_devs[t][k] ~ student_t(nu[1], 0, 1);
      }
    }
  }

  // Likelihood
  if (family == 1) {
    for (i in 1:N) y[i] ~ normal(eta[y_indx[i]], phi[1]);
  }
  if (family == 2) {
    for (i in 1:N) y_int[i] ~ bernoulli_logit(eta[y_indx[i]]);
  }
  if (family == 3) {
    for (i in 1:N) y_int[i] ~ poisson_log(eta[y_indx[i]]);
  }
  if (family == 4) {
    for (i in 1:N) y_int[i] ~ neg_binomial_2_log(eta[y_indx[i]], phi[1]);
  }
  if (family == 5) {
    for (i in 1:N) y[i] ~ gamma(phi[1], phi[1] ./ exp(eta[y_indx[i]]));
  }
  if (family == 6) {
    for (i in 1:N) y[i] ~ lognormal(eta[y_indx[i]], phi[1]);
  }
}
generated quantities {
  vector[N] log_lik;
  if(family==1) {
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | eta[y_indx[n]], phi[1]); // Gaussian
  }
  if(family==2) {
    for (n in 1:N) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(eta[y_indx[n]])); // binomial
  }
  if(family==3) {
    for (n in 1:N) log_lik[n] = poisson_lpmf(y_int[n] | exp(eta[y_indx[n]])); // Poisson
  }
  if(family==4) {
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | eta[y_indx[n]], phi[1]); // NegBin2
  }
  if(family==5) {
    for (n in 1:N) log_lik[n] = gamma_lpdf(y[n] | phi[1], phi[1] ./ exp(eta[y_indx[n]])); // Gamma
  }
  if(family==6) {
    for (n in 1:N) log_lik[n] = lognormal_lpdf(y[n] | eta[y_indx[n]], phi[1]); // Lognormal
  }
}
