# Script for tau-leaping stochastic model 

library(adaptivetau)
library(dplyr)
library(deSolve)
library(ggplot2)
library(patchwork)

### Transitions & Rates ----

transitions = list(c(S = -1, E = +1),# exposure
                  c(E = -1, I = +1), # infection
                  c(I = -1, R = +1)) # recovery
                

lvrates <- function(x, params, t){
  return(c(
    params$beta* (1 - params$V) * x["I"] /params$N * x["S"], #Rate for event 1
    params$sigma * x["E"], # Rate for event 2
    params$gamma * x["I"] # Rate for event 3
  ))
}

### Base parameter ----
sigma_val <- 1/2
gamma_val <- 1/5
R0_val <- 6

beta_val  <- R0_val * gamma_val

params <- list(
  sigma = sigma_val,
  gamma = gamma_val,
  R0 = R0_val,
  N = 2500,
  V = 0,
  beta = beta_val
)

### Vaccine Efficacy loop ----
vac_eff <- c(0, 0.6, 0.7, 0.8, 0.9, 0.95)

n_iter <- 1000

results_list <- list()

counter <- 1

set.seed(626)
# running multiple iterations:
for( v in vac_eff) {
  
  params$V <- v
  
  for (i in 1:n_iter) {
    
    r = ssa.adaptivetau(
      c(S=2499, E=0, I=1, R=0), # initial pop counts
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
    temp_df$R0 <- R0_val
    temp_df$sigma <- sigma_val
    temp_df$gamma <- gamma_val
    
    # Store in list
    results_list[[counter]] <- temp_df
    counter <- counter + 1
  }
}

# Combine all list elements into one large data frame
all_sim_data <- bind_rows(results_list)

LTP_sim_data <- all_sim_data %>% filter(time == 200)

#write.csv(LTP_sim_data, file = "Baseline_model_LTP.csv", row.names = F)
saveRDS(all_sim_data, file = "Baseline_model.rds")

### Averaging tau-leaping sims

time_grid <- seq(from= 0, to = 200, by = 0.5)

interpolated_data <- all_sim_data %>%
  group_by(vaccine_efficacy, run_number) %>%
  reframe(
    
    S = approx(x = time, y = S, xout = time_grid, method = "linear", rule = 2)$y,
    E = approx(x = time, y = E, xout = time_grid, method = "linear", rule = 2)$y,
    I = approx(x = time, y = I, xout = time_grid, method = "linear", rule = 2)$y,
    R = approx(x = time, y = R, xout = time_grid, method = "linear", rule = 2)$y,
    
    time = time_grid
  )

prev_results <- interpolated_data %>%
  group_by(vaccine_efficacy, time) %>%
  summarise(
    
    #Median 
    median_I = median(I),
    u_I = quantile(I, probs = 0.975),
    l_I = quantile(I, probs = 0.025)
    
    ) %>%
  mutate( 
    I = median_I / 2500 * 100,
    u_I = u_I / 2500 * 100,
    l_I = l_I / 2500 * 100)
  

plot_list <- list() #storage of individual plots

efficacies <- unique(prev_results$vaccine_efficacy)

# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- prev_results %>% filter(vaccine_efficacy == v)
  
  # plot generation
  p <- ggplot(data_subset, aes(x = time, y = I))+
    geom_ribbon(aes(ymin = l_I, ymax = u_I), alpha = 0.2, color = "steelblue", fill = "steelblue")+
    geom_line(linewidth = 1)+
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = "Time (Days)",
         y = "Prevalence (%)")+
    coord_cartesian(ylim = c(0, max(prev_results$u_I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3)+
  plot_annotation(
    title = 'Stochastic Infection Trajectories by Vaccine Efficacy')

print(final_plot)

ggsave("Figure_1.tiff", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)


# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- averaged_results %>% filter(vaccine_efficacy == v)
  
  # plot generation
  p <- ggplot(data_subset, aes(x = time, y = mean_I))+
    geom_line(color = "black", linewidth = 1.2)+
    geom_line(data = data_subset, aes(x = time, y = l_I), color = "maroon", linetype = "dashed", linewidth = 1)+
    geom_line(data = data_subset, aes(x = time, y = u_I), color = "maroon", linetype = "dashed", linewidth = 1)+
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = "Time (Days)",
         y = "Prevalence (%)")+
    coord_cartesian(ylim = c(0, max(averaged_results$I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3)+
  plot_annotation(
    title = 'Stochastic Infection Trajectories by Vaccine Efficacy')

print(final_plot)

ggsave("test_run_mean.png", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)


### Attack rate ----
N= 2500 # size of the farm

cum_df <- all_sim_data %>%
  filter(time == 200) %>%
  mutate(Cum_Inf = N - S) 
  
cum_stats <- cum_df %>%
  group_by(vaccine_efficacy, time) %>%
  summarise(
    
    # Mean
    Mean_Cum = mean(Cum_Inf),
    
    # Quantile
    quantile_lower = quantile(Cum_Inf, 0.025),
    quantile_upper = quantile(Cum_Inf, 0.975),
    
    # SE
    sd_cum = sd(Cum_Inf),
    se_cum = sd_cum / sqrt(2500),
    
    # Normal 95% CI
    ci_lower = Mean_Cum - (1.96 * se_cum),
    ci_upper = Mean_Cum + (1.96 * se_cum),
    
    # Poisson 95% CI
    p_ci_lower = Mean_Cum - (1.96 * sqrt(Mean_Cum/2500)),
    p_ci_upper = Mean_Cum + (1.96 * sqrt(Mean_Cum/2500)),
    
    .groups = "drop"
    )

AR_stats <- cum_df %>%
  mutate(AR = Cum_Inf / N * 100) %>%
  group_by(vaccine_efficacy, time) %>%
  summarise(
    
    # Mean
    Mean_AR = mean(AR),
    
    # SE
    sd_AR = sd(AR),
    se_AR = sd_AR / sqrt(N),
    
    # Normal 95% CI
    ci_lower = Mean_AR - (1.96 * se_AR),
    ci_upper = Mean_AR + (1.96 * se_AR),
    
    .groups = "drop"
  )


hist_list <- list()
 efficacies <- unique(cum_df$vaccine_efficacy)
 
for ( v in efficacies){
  
   hist_subset <- cum_df %>% filter(vaccine_efficacy == v)
   
   p <- ggplot(data = hist_subset, aes(x= Cum_Inf))+
     geom_histogram(bins = 100)
   
   hist_list[[as.character(v)]] <- p
 }
 
 hist_plot <- wrap_plots(hist_list, ncol = 3)
 
 print(hist_plot)
### VISUALIZATION Mean ----

plot_list <- list() #storage of individual plots

efficacies <- unique(all_sim_data$vaccine_efficacy)

# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- all_sim_data %>% filter(vaccine_efficacy == v)
  
  # # filter mean efficacy
  # mean_subset <- averaged_results %>% filter(vaccine_efficacy == v)
  
  # plot generation
  p <- ggplot()+
    geom_line(data = data_subset,
              aes(x = time, y = I, group=run_number, color = "Individual Run"),alpha = 0.1)+
    scale_color_manual(name = NULL, # Set to "Legend" if you want a title
                       values = c("Individual Run" = "steelblue")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = "Time (Days)",
         y = "Number of Cases")+
    coord_cartesian(ylim = c(0, max(all_sim_data$I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3) +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom")
#   plot_annotation(
#     title = 'Stochastic Infection Trajectories by Vaccine Efficacy')

print(final_plot)

ggsave("Epi_curves.tiff", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)

### Deterministic model ----
# Model inputs
N <- 2500  # herd size
sigma <- 1/2 # latent 
gamma <- 1/5  # infectious rate

time <- seq(0, 200, by = 1)
initial_state <- c(S=2499, E=0, I=1, R=0)
R0 <- 6

determinstic_result_list <- list()

# SEIR model function
seir_model <- function(time, state, parameters) {
  S <- state["S"]
  E <- state["E"]
  I <- state["I"]
  R <- state["R"]
  
  R0<- parameters["R0"]
  V <- parameters["V"]
  sigma <- parameters["sigma"]
  gamma <- parameters["gamma"]
  mu <- parameters["mu"]
  N <- parameters["N"]
  
  beta <- R0*gamma
  
  dS <- -beta* (1-V) * S * I /N
  dE <- beta * (1-V) * S * I / N - sigma * E
  dI <- sigma * E- gamma * I - mu * I
  dR <- gamma * I

  list(c(dS, dE, dI, dR))
}

for (v in vac_eff) {
  
  parameters <- c(
    sigma = sigma,
    gamma = gamma,
    R0 = R0,
    N = N,
    V = v
  )
  
  # Solve SEIR for vaccinated
  out <- ode(y = initial_state, times = time, func = seir_model, parms = parameters,
             rtol = 1e-6, atol = 1e-8)
  
  temp_df <- as.data.frame(out)
  
  temp_df$vaccine_efficacy <- v
  
  determinstic_result_list[[as.character(v)]] <- temp_df
}

all_deterministic_data <- bind_rows(determinstic_result_list)

write.csv(all_deterministic_data, file = "Baseline_deterministic_model.csv", row.names = F)

### VISUALIZATION  Deterministic ----
plot_list <- list() #storage of individual plots

efficacies <- unique(all_sim_data$vaccine_efficacy)

# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- all_sim_data %>% filter(vaccine_efficacy == v)
  
  # filter mean efficacy
  det_subset <- all_deterministic_data %>% filter(vaccine_efficacy == v)
  
  # plot generation
  p <- ggplot(data_subset, aes(x = time, y = I))+
    geom_line(aes(group=run_number), color = "steelblue", alpha = 0.1)+
    geom_line(data = det_subset, aes(x = time, y = I), color = "black", linewidth = 1)+
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = NULL,
         y = NULL)+
    coord_cartesian(ylim = c(0, max(all_sim_data$I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3)+
  plot_annotation(
    title = 'Stochastic Infection Trajectories by Vaccine Efficacy')

print(final_plot)

ggsave("test_run.png", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)



