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
}
transformed data {
  vector[n_varying_covars] zeros;
  for(i in 1:n_varying_covars) zeros[i] = 0;
}
parameters {
  vector[n_fixed_covars] b_fixed;
  vector[n_varying_covars] b_devs0;
  vector[n_varying_covars] b_devs[nT-1];
  cholesky_factor_corr[n_varying_covars * correlated_rw] Lcorr; // cholesky factor (L_u matrix for R)
  vector<lower=0>[n_varying_covars] sigma;
  vector<lower=0>[1] phi;
  vector<lower=0>[est_df] nu;
  vector[n_fixed_NAs] missing_fixed;
  vector[n_varying_NAs] missing_varying;
}
transformed parameters {
  vector[n_varying_covars] b_varying[nT];
  //corr_matrix[n_varying_covars * correlated_rw] R; // correlation matrix
  //cov_matrix[n_varying_covars * correlated_rw] Sigma; // VCV matrix derived
  matrix[n_varying_covars * correlated_rw, n_varying_covars * correlated_rw] R;
  matrix[n_varying_covars * correlated_rw, n_varying_covars * correlated_rw] Sigma;
  vector[nT] eta; // linear predictor on link scale
  matrix[nT, n_fixed_covars] X_fixed;
  matrix[nT, n_varying_covars] X_varying;

  // covariance matrix stuff
  if(correlated_rw == 1) {
    R = multiply_lower_tri_self_transpose(Lcorr); // R = Lcorr * Lcorr'
    Sigma = quad_form_diag(R, sigma); // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
  }

  // re-assemble X matrices for fixed and time-varying effects
  // for(i in 1:nT) {
  //   for(j in 1:n_fixed_covars) {
  //     X_fixed[i,j] = 0;
  //   }
  //   for(j in 1:n_varying_covars) {
  //     X_varying[i,j] = 0;
  //   }
  // }
  if(n_fixed_covars > 0) {
    for(i in 1:fixed_N) {
      X_fixed[fixed_time_indx[i], fixed_var_indx[i]] = fixed_x_value[i];
    }
    // add in missing vals
    for(i in 1:n_fixed_NAs) {
      X_fixed[fixed_time_indx[fixed_NAs[i]], fixed_var_indx[fixed_NAs[i]]] = missing_fixed[i];
    }
  }
  if(n_varying_covars > 0) {
    for(i in 1:varying_N) {
      X_varying[varying_time_indx[i], varying_var_indx[i]] = varying_x_value[i];
    }
    // add in missing vals
    for(i in 1:n_varying_NAs) {
      X_varying[varying_time_indx[varying_NAs[i]], varying_var_indx[varying_NAs[i]]] = missing_varying[i];
    }
  }

  // time - varying coefficients
  b_varying[1] = b_devs0;
  for(t in 2:nT) {
    b_varying[t] = b_varying[t-1] + b_devs[t-1];
  }

  // calculate predictions (eta)
  for(t in 1:nT) eta[t] = 0;
  if(n_fixed_covars > 0) eta = eta + X_fixed * b_fixed;
  if(n_varying_covars > 0) {
    for(t in 1:nT) {
      eta[t] = eta[t] + X_varying[t] * b_varying[t];
    }
  }
}
model {
  sigma ~ cauchy(0, 5); // prior for sigma
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  phi ~ student_t(3,0,2); // obseervation variance
  b_fixed ~ normal(0,1);
  nu ~ student_t(3,0,2);

  missing_fixed ~ normal(0,1); // estimates of missing Xs for fixed model
  missing_varying ~ normal(0,1); // estimates of missing Xs for time varying model
  b_devs0 ~ normal(0,1); // initial values of B at time t
  if(est_df == 0) {
    if(correlated_rw == 1) {
      for(t in 1:(nT-1)) {
        b_devs[t] ~ multi_normal(zeros, Sigma);
      }
    } else {
      for(t in 1:(nT-1)) {
        b_devs[t] ~ normal(zeros, sigma);
      }
    }
  } else {
    if(correlated_rw == 1) {
      for(t in 1:(nT-1)) {
        b_devs[t] ~ multi_student_t(nu[1], zeros, Sigma);
      }
    } else {
      for(t in 1:(nT-1)) {
        b_devs[t] ~ student_t(nu[1], zeros, sigma);
      }
    }
  }

  if(family==1) {
    for (i in 1:N) y[i] ~ normal(eta[y_indx[i]], phi[1]);
    //y ~ normal(0,1);//normal(eta, phi[1]); // Gaussian
  }
  if(family==2) {
    for (i in 1:N) y_int[i] ~ bernoulli_logit(eta[y_indx[i]]);
    //y_int ~ bernoulli_logit(eta); // binomial
  }
  if(family==3) {
    for (i in 1:N) y_int[i] ~ poisson_log(eta[y_indx[i]]);
    //y_int ~ poisson_log(eta); // Poisson
  }
  if(family==4) {
    for (i in 1:N) y_int[i] ~ neg_binomial_2_log(eta[y_indx[i]], phi[1]);
    //y_int ~ neg_binomial_2_log(eta, phi[1]); // NegBin2
  }
  if(family==5) {
    for (i in 1:N) y[i] ~ gamma(phi[1], phi[1] ./ exp(eta[y_indx[i]]));
    //y ~ gamma(phi[1], phi[1] ./ exp(eta)); // Gamma
  }
  if(family==6) {
    for (i in 1:N) y[i] ~ lognormal(eta[y_indx[i]], phi[1]);
    //y ~ lognormal(eta, phi[1]); // Lognormal
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
