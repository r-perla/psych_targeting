// Meta-analysis model combining machine learning prediction accuracy
// with advertising targeting effects across personality dimensions

data {
  // Effect and study counts
  int n_ad_effects;          // Number of advertising targeting effect sizes
  int n_ml_effects;          // Total number of machine learning effect sizes
  int n_ml_studies;          // Number of unique ML studies for random effects
  
  // Raw data vectors
  vector[n_ad_effects] ad_effects;  // Standardized effects from ad targeting
  vector[n_ad_effects] ad_ses;      // Standard errors of ad effects
  vector[n_ml_effects] ml_r_squared;  // R² from ML prediction studies
  vector[n_ml_effects] ml_ns;         // Sample sizes for ML studies
  
  // Grouping variables
  array[n_ml_effects] int ml_study_id;  // Study ID for random effects (maps effects to studies)
  array[n_ml_effects] int dimension;     // Personality dimension (1-5) for each effect
}

transformed data {
  // Convert ML metrics to Fisher's z scale for meta-analysis
  vector[n_ml_effects] ml_r = sqrt(ml_r_squared);  // Convert R² to r
  vector[n_ml_effects] ml_z = 0.5 * log((1 + ml_r) ./ (1 - ml_r));  // r to z
  
  // Convert ad effects to Fisher's z scale
  vector[n_ad_effects] ad_z = .5 * log((1 + ad_effects) ./ (1 - ad_effects));
  
  // Calculate standard errors for ML effects using standard formula
  vector[n_ml_effects] ml_ses = 1 / sqrt(ml_ns - 3);
}

parameters {
  // Population-level parameters
  real mu_ml;              // Mean Fisher's z for ML studies
  real mu_ad;              // Mean Fisher's z for ad studies
  
  // Variance components
  real<lower=0> tau_ml;    // Between-study SD for ML
  real<lower=0> tau_ad;    // Between-study SD for ads
  
  // Degrees of freedom for robustness
  real<lower=1> nu_ad;     // df for ad studies t-distribution
  real<lower=1> nu_ml;     // df for ML studies t-distribution
  
  // Random and fixed effects
  vector[n_ml_studies] u_ml;  // Study-level random effects
  vector[4] personality_effects;   // Fixed effects for dimensions 2-5
  
  // Effect-specific deviations
  vector[n_ml_effects] theta_ml;   // ML effect-specific deviations
  vector[n_ad_effects] theta_ad;   // Ad effect-specific deviations
}

model {
  // Prior distributions
  mu_ml ~ normal(0.2, 0.1);           // Expect positive ML accuracy
  mu_ad ~ normal(0, 0.5);             // Neutral prior on ad effects
  tau_ml ~ cauchy(0, 5);              // Conservative heterogeneity prior
  tau_ad ~ cauchy(0, 5);              // Conservative heterogeneity prior
  nu_ad ~ gamma(2, 0.1);              // df prior for robustness
  nu_ml ~ gamma(2, 0.1);              // df prior for robustness
  u_ml ~ normal(0, 0.1);              // Random effects prior
  personality_effects ~ normal(0, 0.1); // Personality effects prior
  
  // Construct full personality effects vector (including baseline)
  vector[5] full_personality_effects;
  full_personality_effects[1] = 0;  // baseline dimension
  full_personality_effects[2:5] = personality_effects;
  
  // Hierarchical model for ML studies
  theta_ml ~ student_t(nu_ml, 
                      mu_ml + u_ml[ml_study_id] + full_personality_effects[dimension], 
                      tau_ml);
  ml_z ~ student_t(nu_ml, theta_ml, ml_ses);
  
  // Hierarchical model for ad studies
  theta_ad ~ student_t(nu_ad, mu_ad, tau_ad);
  ad_z ~ student_t(nu_ad, theta_ad, ad_ses);
}

generated quantities {
  // Convert results back to correlation scale
  real mu_r_ml = (exp(2*mu_ml) - 1) / (exp(2*mu_ml) + 1);
  real mu_ad_r = (exp(2*mu_ad) - 1) / (exp(2*mu_ad) + 1);
  
  // Calculate personality-specific correlations
  vector[5] personality_r;
  personality_r[1] = (exp(2*(mu_ml)) - 1) / (exp(2*(mu_ml)) + 1);
  for (d in 2:5) {
    personality_r[d] = (exp(2*(mu_ml + personality_effects[d-1])) - 1) / 
                      (exp(2*(mu_ml + personality_effects[d-1])) + 1);
  }
  
  // Store full personality effects for output
  vector[5] full_personality_effects;
  full_personality_effects[1] = 0;
  full_personality_effects[2:5] = personality_effects;
  
  // Calculate summary effects
  real total_mu_r = mean(personality_r);
  real attenuated_effect = total_mu_r * mu_ad_r;
  
  // Calculate I² heterogeneity statistics
  real i_squared_ad = square(tau_ad) / (square(tau_ad) + mean(square(ad_ses))) * 100;
  real i_squared_ml = square(tau_ml) / (square(tau_ml) + mean(square(ml_ses))) * 100;
  
  // Compute final estimates combining fixed and random effects
  vector[n_ml_studies] final_theta_z = theta_ml + u_ml[ml_study_id] + 
                                      full_personality_effects[dimension];
  
  // Convert final estimates to correlation scale
  vector[n_ml_studies] final_theta_r;
  for (i in 1:n_ml_studies) {
    final_theta_r[i] = (exp(2*final_theta_z[i]) - 1) / 
                       (exp(2*final_theta_z[i]) + 1);
  }
}
