data {
  int n_ad_studies; // Number of targeting studies
  int n_ml_studies; // Number of ML studies
  int n_ind_ml_studies; // Unique ML studies
  vector[n_ad_studies] ad_effects;
  vector[n_ad_studies] ad_ses;
  vector[n_ml_studies] ml_r_squared;
  vector[n_ml_studies] ml_ns;
  array[n_ml_studies] int ml_study_id;
  array[n_ml_studies] int dimension;
}

transformed data {
  vector[n_ml_studies] ml_r  = sqrt(ml_r_squared);
  vector[n_ml_studies] ml_z  = 0.5 * log((1 + ml_r) ./ (1 - ml_r));
  vector[n_ad_studies] ad_z = 0.5 * log((1 + ad_effects) ./ (1 - ad_effects));
  vector[n_ml_studies] ml_ses = 1 ./ sqrt(ml_ns - 3);
}

parameters {
  real mu_ml;
  real mu_ad;
  real<lower=0> tau_ml;
  real<lower=0> tau_ad;
  real<lower=0> sigma_u;

  vector[4] personality_effects;

  // raw parameters for non-centered reparameterization
  vector[n_ind_ml_studies] u_ml_raw;
  vector[n_ml_studies]     theta_ml_raw;
  vector[n_ad_studies]     theta_ad_raw;
}

transformed parameters {
  // transform raw to actual random effects
  vector[n_ind_ml_studies] u_ml = sigma_u * u_ml_raw;
  vector[n_ml_studies] theta_ml;
  vector[n_ad_studies] theta_ad;

  // build full personality vector
  vector[5] full_personality_effects;
  full_personality_effects[1] = 0;
  full_personality_effects[2:5] = personality_effects;

  for (i in 1:n_ml_studies) {
    theta_ml[i] = mu_ml
                + u_ml[ml_study_id[i]]
                + full_personality_effects[dimension[i]]
                + tau_ml * theta_ml_raw[i];
  }
  
  for (j in 1:n_ad_studies) {
    theta_ad[j] = mu_ad + tau_ad * theta_ad_raw[j];
  }
}

model {
  // Priors on population-level
  mu_ml ~ normal(0.2, 0.5);
  mu_ad ~ normal(0, 0.5);
  tau_ml ~ student_t(4, 0, 0.3);
  tau_ad ~ student_t(4, 0, 0.3);
  sigma_u ~ normal(0, 0.3);
  personality_effects ~ normal(0, 0.1);

  // Priors on raw deviations
  u_ml_raw ~ normal(0, 1);
  theta_ml_raw ~ normal(0, 1);
  theta_ad_raw ~ normal(0, 1);

  // Likelihoods
  ml_z ~ normal(theta_ml, ml_ses);
  ad_z ~ normal(theta_ad, ad_ses);
}

generated quantities {
  real mu_r_ml = (exp(2*mu_ml) - 1) / (exp(2*mu_ml) + 1);
  real mu_ad_r = (exp(2*mu_ad) - 1) / (exp(2*mu_ad) + 1);

  vector[5] personality_r;
  personality_r[1] = mu_r_ml;
  for (d in 2:5) {
    real tmp = mu_ml + personality_effects[d-1];
    personality_r[d] = (exp(2*tmp) - 1) / (exp(2*tmp) + 1);
  }

  real total_mu_r = tanh(mu_ml);
  real attenuated_effect    = total_mu_r * mu_ad_r;
  real disattenuated_effect = (abs(total_mu_r) < 1e-6) ? 0 : mu_ad_r / total_mu_r;

  real i_squared_ad = square(tau_ad) / (square(tau_ad) + mean(square(ad_ses))) * 100;
  real i_squared_ml = square(tau_ml) / (square(tau_ml) + mean(square(ml_ses))) * 100;

  vector[n_ml_studies] final_theta_r = tanh(theta_ml);
}
