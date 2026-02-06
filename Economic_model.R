# Economic Model

library(dplyr)
library(progress)

# Load Tau-leap results
SEIR_simulation <- readRDS("Baseline_model.rds")

# Filter the data for the final time point & count number of sick pigs
ftp_data <- SEIR_simulation %>% 
  filter(time==200) %>%
  mutate(Cum_inf = I + R)

# Separate Baseline and Intervention
baseline_data <- ftp_data %>% filter(vaccine_efficacy==0)

intervention_data <- ftp_data %>% filter(vaccine_efficacy > 0)

# Mean baseline results
baseline_mean_inf <- mean(baseline_data$Cum_inf)

# Inputs for vaccination parameters
vaccine_strat_set_cost <- c(1.34, 1.68, 2.02, 2.36, 2.70)
N <- 2500


# Initializing run
econ_results_list <- data.frame()

set.seed(626)
n_sim <- 100

#Setting up a progress bar
total_runs <- length(vaccine_strat_set_cost) * nrow(intervention_data) * n_sim
pb <- progress_bar$new( format = "Running [:bar] :percent ETA: :eta",
                        total= total_runs, clear = F, width = 60)

### Economic Sensitivity Analysis (OFAT) ----

# 1. Define Parameter Estimates (Baseline, Low, High)
# We store these in a list so we can access them programmatically
econ_params_info <- list(
  penalty_i        = list(base = 4.93,   low = 2.3,    high = 9.8),
  feed_cost_per_kg = list(base = 0.29,   low = 0.19,   high = 0.39),
  mortality        = list(base = 0.022,  low = 0.015,  high = 0.033),
  raising_cost     = list(base = 174.36, low = 152.24, high = 203.50),
  profit_loss      = list(base = 12.86,  low = 3.35,   high = 29.26),
  carcass_disposal = list(base = 9.07,   low = 7.34,   high = 10.66)
)

# 2. Build the List of Scenarios
# Start with the Baseline scenario
scenarios <- list(
  list(name = "Baseline", param = "None", type = "Base")
)

# automatically add Low and High scenarios for each parameter
for (param_name in names(econ_params_info)) {
  # Add Low Scenario
  scenarios[[length(scenarios) + 1]] <- list(
    name = paste(param_name, "Low"), 
    param = param_name, 
    type = "Low"
  )
  # Add High Scenario
  scenarios[[length(scenarios) + 1]] <- list(
    name = paste(param_name, "High"), 
    param = param_name, 
    type = "High"
  )
}

# Initialize list to store results (faster than rbind)
results_list <- list()
counter <- 1
pb <- txtProgressBar(min = 0, max = length(vaccine_strat_set_cost) * nrow(intervention_data) * length(scenarios), style = 3)

# 3. Main Loop
for (cost_vax in vaccine_strat_set_cost) {
  
  for (run_number in 1:nrow(intervention_data)) {
    
    # Extract disease data for this specific run
    current_run <- intervention_data[run_number,]
    eff_val <- current_run$vaccine_efficacy
    Total_I <- current_run$Cum_inf
    
    # Loop through the defined sensitivity scenarios (Baseline, Highs, Lows)
    for (scenario in scenarios) {
      
      # A. SET PARAMETER VALUES
      # First, set everything to Baseline
      penalty_i        <- econ_params_info$penalty_i$base
      feed_cost_per_kg <- econ_params_info$feed_cost_per_kg$base
      mortality        <- econ_params_info$mortality$base
      raising_cost     <- econ_params_info$raising_cost$base
      profit_loss      <- econ_params_info$profit_loss$base
      carcass_disposal <- econ_params_info$carcass_disposal$base
      
      # Second, overwrite ONLY the specific parameter being tested (if not baseline)
      if (scenario$type == "Low") {
        assign(scenario$param, econ_params_info[[scenario$param]]$low)
      } else if (scenario$type == "High") {
        assign(scenario$param, econ_params_info[[scenario$param]]$high)
      }
      
      # B. CALCULATE INTERVENTION COSTS
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
      
      
      # C. BASELINE / NO VACCINATION COST
      # (Note: Use the SAME economic parameters to compare fairly)
      
      # Feed Penalty
      baseline_extra_feed_kg <- baseline_mean_inf * penalty_i
      baseline_extra_feed_cost <- baseline_extra_feed_kg * feed_cost_per_kg
      
      # Mortality
      baseline_death <- baseline_mean_inf * mortality
      baseline_mortality_cost <- baseline_death * cost_per_death_unit
      
      baseline_total_cost <- baseline_extra_feed_cost + baseline_mortality_cost
      
      
      # D. ECONOMIC MEASURES
      # Benefit
      feed_benefit = baseline_extra_feed_cost - extra_feed_cost
      mortality_benefit = baseline_mortality_cost - mortality_cost
      
      total_benefit = feed_benefit + mortality_benefit
      
      # BCR
      # Handle division by zero if vax cost is 0
      if(total_vax_cost > 0) {
        BCR = total_benefit / total_vax_cost
      } else {
        BCR = NA 
      }
      
      # NPV
      NPV = total_benefit - total_vax_cost
      
      # ROI
      if(total_vax_cost > 0) {
        ROI = (NPV / total_vax_cost) * 100
      } else {
        ROI = NA
      }
      
      # E. STORE RESULTS
      results_list[[counter]] <- data.frame(
        disease_id = run_number,
        scenario_name = scenario$name, # e.g. "Feed Cost High"
        scenario_type = scenario$type, # e.g. "High"
        changed_param = scenario$param, # e.g. "feed_cost_per_kg"
        
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