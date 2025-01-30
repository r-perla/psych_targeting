rm(list = objects())
library(dplyr)
library(ggplot2)
library(plotly)
library(gridExtra)
library(foreach)
library(doParallel)
library(readr)
set.seed(123)

# Define data generating function
gen_data <- function(n, b0, b1, b2, b3) {
  # Generate trait values and conditions
  trait <- rnorm(n) # Measured without error for simplicity
  condition <- rbinom(n, 1, 0.5) # Assuming participants are assigned at random
  
  # Generate outcome from main effects and interaction effects
  y <- b0 + b1 * trait + b2 * condition + b3 * trait * condition + rnorm(n)
  
  # Return as data frame
  return(data.frame(trait, condition, y))
}

# Define function to run all designs
get_results <- function(data) {
  # Run design one (Trait-Treatment Interaction)
  model <- lm(y ~ trait * condition, data = data)
  
  # Extract coefficient for interaction term
  b3_design1 <- coef(model)["trait:condition"]
  
  # Run design two (Selected Population Design)
  d_design2 <- data[data$trait > quantile(data$trait, 0.75), ] # Preselect audience for high trait values
  
  # Regress outcome on condition
  model <- lm(y ~ condition, data = d_design2)
  
  # Extract coefficient for condition
  b3_design2 <- coef(model)["condition"]
  
  # Run design 3 (Matched-mismatched design)
  data <- data %>%
    mutate(matched = case_when(
      (trait > quantile(trait, .75) & condition == 1) | (trait <= quantile(trait, .25) & condition == 0) ~ 1,
      (trait > quantile(trait, .75) & condition == 0) | (trait <= quantile(trait, .25) & condition == 1) ~ 0,
      TRUE ~ NA_real_ # This should never happen
    ))
  
  # Regress outcome in matching condition
  model <- lm(y ~ matched, data = data)
  
  # Extract coefficient for matched condition
  b3_design3 <- coef(model)["matched"]
  
  # Run design 4 (Trait-only design)
  d_design4 <- data %>%
    filter(condition == 1) # Only present tailored at
  
  # Run regression
  model <- lm(y ~ trait, data = d_design4)
  
  # Extract coefficient for trait
  b3_design4 <- coef(model)["trait"]
  
  # Return all coefficients as a dataframe
  return(data.frame(
    design = c("Trait-Treatment Interaction", "Selected Population Design", "Matched-mismatched Design", "Trait-only Design"),
    b3 = c(b3_design1, b3_design2, b3_design3, b3_design4)
  ))
}

# Define grid search function to search through different combinations of main effects
grid_search <- function(grid, n, b3, n_designs) {
  results <- foreach(i = 1:nrow(grid), .combine = rbind, .packages = c("dplyr")) %dopar% {
    # Initialize results vector for this iteration
    result_row <- numeric(n_designs + ncol(grid))
    
    # Extract main effects from grid
    b0 <- grid[i, "b0"]
    b1 <- grid[i, "b1"]
    b2 <- grid[i, "b2"]
    
    # Generate data
    data <- gen_data(n, b0, b1, b2, b3)
    
    # Run all designs
    design_results <- get_results(data)
    
    # Store results in vector
    result_row[1:n_designs] <- as.numeric(design_results$b3)
    result_row[(n_designs + 1):length(result_row)] <- as.numeric(grid[i, ])
    
    result_row
  }
  
  # Convert to dataframe and add column names
  results <- as.data.frame(results)
  colnames(results) <- c(paste0("b3_design", 1:n_designs), colnames(grid))
  
  return(results)
}

# Set global parameters to run the simulation
n <- 10000000 # Asymptotic sample size (ensure this is big enough so that even subset-designs have asymptotic n)
b3 <- .5
n_designs <- 4
grid <- expand.grid(
  #b0 = seq(-1, 1, by = 0.1),
  b0 = 0,
  b1 = c(-.25, 0, .25),
  b2 = c(-.25, 0, .25)
)

# Register parallel backend
cl <- makeCluster(10) # Do not use all your cores, this get's *very* RAM intensive
clusterExport(cl, c("gen_data", "get_results"))
registerDoParallel(cl)

# Run the grid search
results <- grid_search(grid, n, b3, n_designs)

# Stop backend
stopCluster(cl)

# Print results
head(results)

# Write to file
write_rds(results, "data/simulation_results.rds", compress = "gz")

plot_3d_bias <- function(results, true_b3) {
  # Create matrices for each design
  b1_unique <- unique(results$b1)
  b2_unique <- unique(results$b2)
  
  # Create empty matrices for each design
  z_matrices <- list()
  for(i in 1:4) {
    z_matrices[[i]] <- matrix(NA, nrow = length(b1_unique), ncol = length(b2_unique))
    design_col <- paste0("b3_design", i)
    
    # Fill matrices with bias (estimated - true)
    for(j in seq_along(b1_unique)) {
      for(k in seq_along(b2_unique)) {
        idx <- results$b1 == b1_unique[j] & results$b2 == b2_unique[k]
        z_matrices[[i]][j,k] <- results[idx, design_col] - true_b3
      }
    }
  }
  
  # Create plot
  p <- plot_ly()
  
  # Add surface for each design
  designs <- c("Trait-Treatment Interaction",
               "Selected Population Design",
               "Matched-mismatched Design",
               "Trait-only Design")
  colors <- c("black", "#27687B", "#64277B", "#7B3A27")
  
  for(i in 1:4) {
    p <- add_surface(p,
                     x = b1_unique,
                     y = b2_unique,
                     z = z_matrices[[i]],
                     opacity = ifelse(i == 1, 1, 0.5),
                     colorscale = list(c(0,1), c(colors[i], colors[i])),
                     showscale = FALSE,
                     showlegend = TRUE,
                     name = designs[i],
                     legendgroup = designs[i])
  }
  
  # Add layout details
  p <- layout(p,
              scene = list(
                xaxis = list(title = list(text = "Trait main effect", font = list(size = 20))),
                yaxis = list(title = list(text = "Ad condition main effect", font = list(size = 20))),
                zaxis = list(title = list(text = "Bias in targeting effect", font = list(size = 20)))
              ),
              showlegend = F,
              margin = list(r = 120),  # Add right margin for legend
              title = NULL)
  
  return(p)
}

# Create the plot
fig <- plot_3d_bias(results, b3)
fig
