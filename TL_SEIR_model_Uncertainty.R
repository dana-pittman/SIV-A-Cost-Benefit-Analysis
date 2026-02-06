# Script for tau-leaping stochastic model 

library(adaptivetau)
library(dplyr)

### Transitions & Rates ----

transitions = list(c(S = -1, E = +1),# exposure
                  c(E = -1, I = +1), # infection
                  c(I = -1, R = +1)) # recovery

lvrates <- function(x, params, t){
  return(c(
    params$beta* (1 - params$V) * x["I"] /params$N * x["S"], # Rate for event 1
    params$sigma * x["E"], # Rate for event 2
    params$gamma * x["I"] # Rate for event 3
  ))
}


### SAMPLE SIGMA (LATENT PERIOD)  ----
gamma_val <- 1/5
R0_val <- 6

beta_val <- R0_val * gamma_val

params <- list(
  gamma = gamma_val,
  R0 = R0_val,
  beta = beta_val,
  N = 1000,
  V = 0
)

### Vaccine Efficacy loop ----
vac_eff <- c(0, 0.6, 0.7, 0.8, 0.9, 0.95)

n_iter <- 1000

results_sigma <- list()

counter <- 1

# running multiple iterations:
for( v in vac_eff) {
  
  params$V <- v
  
  for (i in 1:n_iter) {
    
    # Sampling distribution
    current_sigma <-runif(n=1, min = 1/5, max = 1/1.4)
    
    #Update params list
    params$sigma <- current_sigma
    
    #Run Model
    r = ssa.adaptivetau(
      c(S=999, E=0, I=1, R=0), # initial pop counts
      transitions,
      lvrates,
      params,
      tf=200
        )

    # Convert matrix to data frame
    temp_df <- as.data.frame(r)
    
    # Add metadata columns
    temp_df$run_number <- i
    temp_df$vaccine_efficacy <- v
    
    # Analysis info
    temp_df$R0 <- R0_val
    temp_df$sigma <- current_sigma
    temp_df$gamma <- gamma_val
    
    # Store in list
    results_sigma[[counter]] <- temp_df
    counter <- counter + 1
  }
}

# Combine all list elements into one large data frame
df_sigma_uncertainty<- bind_rows(results_sigma)



### SAMPLE GAMMA (INFECTIOUS PERIOD)  ----
sigma_val <- 1/2
R0_val <- 6

params <- list(
  sigma = sigma_val,
  R0 = R0_val,
  N = 1000,
  V = 0,
  beta = beta_val
)

### Vaccine Efficacy loop ---
vac_eff <- c(0, 0.6, 0.7, 0.8, 0.9, 0.95)

n_iter <- 1000

results_gamma <- list()

counter <- 1

# running multiple iterations:
for( v in vac_eff) {
  
  params$V <- v
  
  for (i in 1:n_iter) {
    
    # Sampling distribution
    current_gamma <-runif(n=1, min = 1/10.4, max = 1/2.4)
    
    #Calculating beta
    current_beta <- params$R0 * current_gamma
    
    #Update params list
    params$gamma <- current_gamma
    
    #Run Model
    r = ssa.adaptivetau(
      c(S=999, E=0, I=1, R=0), # initial pop counts
      transitions,
      lvrates,
      params,
      tf=200
    )
    
    # Convert matrix to data frame
    temp_df <- as.data.frame(r)
    
    # Add metadata columns
    temp_df$run_number <- i
    temp_df$vaccine_efficacy <- v
    
    # Analysis info
    temp_df$R0 <- R0_val
    temp_df$sigma <- sigma_val
    temp_df$gamma <- current_gamma
    
    # Store in list
    results_gamma[[counter]] <- temp_df
    counter <- counter + 1
  }
}

# Combine all list elements into one large data frame
df_gamma_uncertainty<- bind_rows(results_gamma)



### SAMPLE R0 ----
sigma_val <- 1/2
gamma_val <- 1/5

params <- list(
  sigma = sigma_val,
  gamma = gamma_val,
  N = 1000,
  V = 0
)

### Vaccine Efficacy loop ----
vac_eff <- c(0, 0.6, 0.7, 0.8, 0.9, 0.95)

n_iter <- 1000

results_R0 <- list()

counter <- 1

# running multiple iterations:
for( v in vac_eff) {
  
  params$V <- v
  
  for (i in 1:n_iter) {
    
    # Sampling distribution
    current_R0 <-runif(n=1, min = 2.5, max = 10.66)
    
    #Calculating beta
    current_beta <- current_R0 * params$gamma
    
    #Update params list
    params$R0 <- current_R0
    params$beta <- current_beta
    
    #Run Model
    r = ssa.adaptivetau(
      c(S=999, E=0, I=1, R=0), # initial pop counts
      transitions,
      lvrates,
      params,
      tf=200
    )
    
    # Convert matrix to data frame
    temp_df <- as.data.frame(r)
    
    # Add metadata columns
    temp_df$run_number <- i
    temp_df$vaccine_efficacy <- v
    
    # Analysis info
    temp_df$R0 <- current_R0
    temp_df$sigma <- sigma_val
    temp_df$gamma <- gamma_val
    
    # Store in list
    results_R0[[counter]] <- temp_df
    counter <- counter + 1
  }
}

# Combine all list elements into one large data frame
df_R0_uncertainty<- bind_rows(results_R0)


uncertainty_all_data <- bind_rows(df_R0_uncertainty, df_gamma_uncertainty, df_sigma_uncertainty)

uncertainty_all_data <- uncertainty_all_data %>% filter(time == 200)

saveRDS(uncertainty_all_data, "Uncertainty_Analysis_alltime.rds")

#write.csv(uncertainty_all_data, "Uncertainty_Analysis.csv", row.names = F)


