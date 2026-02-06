# Script for OWSA using deterministic model

library(ggplot2)
library(dplyr)
library(deSolve)
library(progress)


### Deterministic  Model ----

# Model inputs
N <- 2500  # herd size
time <- seq(0, 200, by = 1)
initial_state <- c(S=2499, E=0, I=1, R=0)

# Define the range of Vaccine Efficacies to test
vac_eff <- c(0, 0.6, 0.7, 0.8, 0.9, 0.95)

# Baseline parameters here
baseline_params <- list(
  sigma = 1/2,    # Baseline latent rate
  gamma = 1/5,    # Baseline infectious rate
  R0 = 6,         # Baseline R0
  N = N,
  mu = 0          # Mortality rate
)

# Scenario changing from the baseline (Epi Parameter Scenarios)
epi_scenarios <- list(
  list(name = "Baseline (R0=6)", param = "R0", value = 6),
  
  # R0 variations
  list(name = "R0 = 2.5",   param = "R0", value = 2.5),
  list(name = "R0 = 10.66", param = "R0", value = 10.66),
  
  # Sigma variations
  list(name = "Sigma = 1/5",   param = "sigma", value = 1/5),
  list(name = "Sigma = 1/1.4", param = "sigma", value = 1/1.4),
  
  # Gamma variations
  list(name = "Gamma = 1/10.4", param = "gamma", value = 1/10.4),
  list(name = "Gamma = 1/2.4",  param = "gamma", value = 1/2.4)
)

deterministic_result_list <- list()

# SEIR model function
seir_model <- function(time, state, parameters) {
  S <- state["S"]
  E <- state["E"]
  I <- state["I"]
  R <- state["R"]
  
  R0    <- parameters["R0"]
  V     <- parameters["V"]
  sigma <- parameters["sigma"]
  gamma <- parameters["gamma"]
  mu    <- parameters["mu"]
  N     <- parameters["N"]
  
  beta <- R0 * gamma
  
  dS <- -beta * (1 - V) * S * I / N
  dE <- beta * (1 - V) * S * I / N - sigma * E
  dI <- sigma * E - gamma * I - mu * I
  dR <- gamma * I
  
  list(c(dS, dE, dI, dR))
}

# Nested Loop: Epi Scenarios x Vaccine Efficacy
counter <- 1 

for (scenario in epi_scenarios) {
  
  # Inner loop: Iterate through vaccine efficacies for THIS scenario
  for (v in vac_eff) {
    
    # Reset to baseline
    current_params <- unlist(baseline_params)
    
    # 1. Apply the epi parameter change
    current_params[scenario$param] <- scenario$value
    
    # 2. Apply the vaccine efficacy
    current_params["V"] <- v
    
    # Solve SEIR
    out <- ode(y = initial_state, 
               times = time, 
               func = seir_model, 
               parms = current_params,
               rtol = 1e-6, atol = 1e-8)
    
    temp_df <- as.data.frame(out)
    
    # Add metadata
    temp_df$scenario_name    <- scenario$name
    temp_df$changed_param    <- scenario$param
    temp_df$param_value      <- scenario$value
    temp_df$vaccine_efficacy <- v
    
    # Store result
    deterministic_result_list[[counter]] <- temp_df
    counter <- counter + 1
  }
}

# Combine all results
all_deterministic_data <- bind_rows(deterministic_result_list)

# Filter for the final timestep to get cumulative results
det_data <- all_deterministic_data %>% 
  filter(time == 200) %>%
  mutate(Cum_inf = I + R) %>%
  select(scenario_name, changed_param, param_value, vaccine_efficacy, Cum_inf)

### Data Bridging ----

# Identify what the "No Vaccination" infection count was for EACH epi scenario (R0=2.5, R0=6, etc.) independently.

# Create a lookup table for the baseline infections (VE = 0)
baseline_lookup <- det_data %>%
  filter(vaccine_efficacy == 0) %>%
  select(scenario_name, baseline_inf_count = Cum_inf)

# Merge this back into the main data
# Now every row knows what the "No Vax" result was for its specific scenario
intervention_data <- det_data %>%
  left_join(baseline_lookup, by = "scenario_name")


### Economic Sensitivity Analysis (OFAT) ----

# Inputs for vaccination parameters
vaccine_strat_set_cost <- c(1.34, 1.68, 2.02, 2.36, 2.70)

# Define Parameter Estimates (Baseline, Low, High)
econ_params_info <- list(
  penalty_i        = list(base = 4.93,   low = 2.3,    high = 9.8),
  feed_cost_per_kg = list(base = 0.29,   low = 0.19,   high = 0.39),
  mortality        = list(base = 0.022,  low = 0.015,  high = 0.033),
  raising_cost     = list(base = 174.36, low = 152.24, high = 203.50),
  profit_loss      = list(base = 12.86,  low = 3.35,   high = 29.26),
  carcass_disposal = list(base = 9.07,   low = 7.34,   high = 10.66)
)

# Build the List of Economic Scenarios
econ_scenarios <- list(
  list(name = "Baseline Econ", param = "None", type = "Base")
)

# automatically add Low and High scenarios for each parameter
for (param_name in names(econ_params_info)) {
  econ_scenarios[[length(econ_scenarios) + 1]] <- list(
    name = paste(param_name, "Low"), 
    param = param_name, 
    type = "Low"
  )
  econ_scenarios[[length(econ_scenarios) + 1]] <- list(
    name = paste(param_name, "High"), 
    param = param_name, 
    type = "High"
  )
}

# Initialize list to store results
results_list <- list()
counter <- 1
total_iterations <- length(vaccine_strat_set_cost) * nrow(intervention_data) * length(econ_scenarios)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

# Main Loop
for (cost_vax in vaccine_strat_set_cost) {
  
  for (run_number in 1:nrow(intervention_data)) {
    
    # Extract disease data for this specific run
    current_run <- intervention_data[run_number,]
    
    eff_val <- current_run$vaccine_efficacy
    Total_I <- current_run$Cum_inf
    
    # CRITICAL: Use the scenario-specific baseline we calculated in Step 2
    baseline_mean_inf <- current_run$baseline_inf_count 
    
    # Loop through the defined sensitivity scenarios (Baseline, Highs, Lows)
    for (scenario in econ_scenarios) {
      
      # SET PARAMETER VALUES
      # First, set everything to Baseline
      penalty_i        <- econ_params_info$penalty_i$base
      feed_cost_per_kg <- econ_params_info$feed_cost_per_kg$base
      mortality        <- econ_params_info$mortality$base
      raising_cost     <- econ_params_info$raising_cost$base
      profit_loss      <- econ_params_info$profit_loss$base
      carcass_disposal <- econ_params_info$carcass_disposal$base
      
      # Overwrite ONLY the specific parameter being tested
      if (scenario$type == "Low") {
        assign(scenario$param, econ_params_info[[scenario$param]]$low)
      } else if (scenario$type == "High") {
        assign(scenario$param, econ_params_info[[scenario$param]]$high)
      }
      
      # CALCULATE INTERVENTION COSTS
      # Feed penalty
      extra_feed_kg <- Total_I * penalty_i
      extra_feed_cost <- extra_feed_kg * feed_cost_per_kg
      
      # Vaccine cost
      total_vax_cost <- N * cost_vax
      
      # Mortality
      num_deaths <- Total_I * mortality
      cost_per_death_unit <- raising_cost + profit_loss + carcass_disposal
      mortality_cost <- num_deaths * cost_per_death_unit
      
      # Total intervention cost
      total_cost_intervention <- extra_feed_cost + total_vax_cost + mortality_cost
      
      
      # BASELINE / NO VACCINATION COST
      # Feed Penalty
      baseline_extra_feed_kg <- baseline_mean_inf * penalty_i
      baseline_extra_feed_cost <- baseline_extra_feed_kg * feed_cost_per_kg
      
      # Mortality
      baseline_death <- baseline_mean_inf * mortality
      baseline_mortality_cost <- baseline_death * cost_per_death_unit
      
      baseline_total_cost <- baseline_extra_feed_cost + baseline_mortality_cost
      
      
      # ECONOMIC MEASURES
      # Benefit
      feed_benefit = baseline_extra_feed_cost - extra_feed_cost
      mortality_benefit = baseline_mortality_cost - mortality_cost
      
      total_benefit = feed_benefit + mortality_benefit
      
      # NPV (Calculated BEFORE ROI)
      NPV = total_benefit - total_vax_cost
      
      # BCR & ROI
      if(total_vax_cost > 0) {
        BCR = total_benefit / total_vax_cost
        ROI = (NPV / total_vax_cost) * 100
      } else {
        BCR = NA 
        ROI = NA
      }
      
      # STORE RESULTS
      results_list[[counter]] <- data.frame(
        # Meta Data from Epi Model
        bio_scenario_name = current_run$scenario_name,
        bio_changed_param = current_run$changed_param,
        
        # Meta Data from Econ Model
        econ_scenario_name = scenario$name, 
        econ_scenario_type = scenario$type, 
        econ_changed_param = scenario$param, 
        
        VE = eff_val,
        vaccine_cost_per_pig = cost_vax,
        
        # Outcomes
        total_intervention_cost = round(total_cost_intervention, 2),
        baseline_cost_avg = round(baseline_total_cost, 2),
        NPV = round(NPV, 2),
        BCR = round(BCR, 2),
        ROI = round(ROI, 2)
      )
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
}
close(pb)

# Combine all results into one dataframe
econ_results_df <- bind_rows(results_list)

# Quick check of the data
head(econ_results_df)