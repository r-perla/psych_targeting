# Closed-form counterparts to the simulation in design_simulations.Rmd
# DGP: y = b0 + b1*T + b2*C + b3*T*C + e, T ~ N(0,1), C in {0,1} randomized, e ~ N(0,1)
# Designs:
#  1) TraitxTreatment (correct estimand): E[hat] = b3
#  2) Selected-population (T > qnorm(p)):  E[hat] = b2 + b3 * lambda_p
#  3) Matched–mismatched (T>qnorm(p) vs T<=-qnorm(p)): E[hat] = b3 * lambda_p
#  4) Trait-only within C=1: E[hat] = b1 + b3
# where lambda_p = phi(qnorm(p)) / (1 - Phi(qnorm(p))) is the upper-tail truncated-normal mean.

# Clear workspace and load required packages
rm(list = objects())
library(dplyr)
library(tidyr)
library(ggplot2)

# Upper-tail truncated-normal mean lambda_p for T ~ N(0,1) and threshold at qnorm(p)
lambda_upper <- function(p) {
  if (!is.numeric(p) || any(p <= 0.5 | p >= 1)) {
    stop("p must be in (0.5, 1); e.g., 0.75 for top quartile.")
  }
  c <- qnorm(p)
  dnorm(c) / (1 - p)
}

# Closed-form expectations for each design
design_estimates <- function(b0, b1, b2, b3, p = 0.75) {
  lam <- lambda_upper(p)
  
  tibble::tibble(
    design = c(
      "Trait×Treatment Interaction",
      "Selected Population",
      "Matched–Mismatched",
      "Trait-only within C=1"
    ),
    estimate = c(
      b3,
      b2 + b3 * lam,
      b3 * lam,
      b1 + b3
    ),
    deviation_from_b3 = estimate - b3
  )
}

# Grid over b1 and b2 (b0 drops out of all expressions)
closed_form_grid <- function(b1_vals, b2_vals, b3, p = 0.75) {
  tidyr::expand_grid(b1 = b1_vals, b2 = b2_vals) %>%
    rowwise() %>%
    do({
      b1 <- .$b1
      b2 <- .$b2
      design_estimates(b0 = 0, b1 = b1, b2 = b2, b3 = b3, p = p) %>%
        mutate(b1 = b1, b2 = b2)
    }) %>%
    ungroup()
}

# Parameters
b3 <- 0.5
p <- 0.75
b1_vals <- c(-0.25, 0, 0.25)
b2_vals <- c(-0.25, 0, 0.25)

# Table of closed-form results
results <- closed_form_grid(b1_vals, b2_vals, b3, p)
print(results)

# Mean absolute deviation from b3 per design (over the b1,b2 grid)
mad_table <- results %>%
  group_by(design) %>%
  summarise(mean_abs_deviation = mean(abs(deviation_from_b3)), .groups = "drop")

cat("\nMean absolute deviation from b3 (over the grid):\n")
print(mad_table)

# 2x2 heatmaps of deviation_from_b3 across (b1,b2)
gg <- results %>%
  mutate(design = factor(
    design,
    levels = c(
      "Trait×Treatment Interaction",
      "Selected Population",
      "Matched–Mismatched",
      "Trait-only within C=1"
    )
  )) %>%
  ggplot(aes(x = b1, y = b2, fill = deviation_from_b3)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", deviation_from_b3)), size = 3) +
  scale_fill_gradient2(
    name = "Deviation\n(estimate - b3)",
    limits = NULL
  ) +
  facet_wrap(~design, ncol = 2) +
  labs(
    x = "Trait main effect (b1)",
    y = "Ad main effect (b2)",
    title = sprintf(
      "Closed-form deviation from personalization effect b3=%.2f (p=%.2f, λ=%.3f)",
      b3, p, lambda_upper(p)
    )
  ) +
  theme_minimal(base_size = 12)

print(gg)

# Note: When b1=0, design 4 coincides with b3; when b2=0, design 2 reduces to b3*lam.
# These equalities are by construction (estimand mismatch, not sampling error).
