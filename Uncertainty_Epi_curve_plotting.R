library(ggplot2)
library(dplyr)
library(patchwork)

#load library
all_data <- readRDS("Uncertainty_Analysis_alltime.rds")

### R0 EPI-CURVE PLOTTING ----

# seperate by varried parameter
R0_varied <- all_data %>% filter(R0 != 6)


# preparing to plot
plot_list <- list() #storage of individual plots

efficacies <- unique(R0_varied$vaccine_efficacy)

# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- R0_varied %>% filter(vaccine_efficacy == v) %>%
    mutate(parm_cat = if_else(R0 > 6, "High", "Low"))

  # plot generation
  p <- ggplot()+
    geom_line(data = data_subset,
              aes(x = time, y = I, group=run_number, color = parm_cat),alpha = 0.1)+
    scale_color_manual(name = NULL, # Set to "Legend" if you want a title
                       values = c("High" = "darkred", "Low" = "steelblue")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = "Time (Days)",
         y = "Number of Cases")+
    coord_cartesian(ylim = c(0, max(R0_varied$I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
   plot_annotation(
     title = 'Stochastic Infection Trajectories by Vaccine Efficacy (R0)')

print(final_plot)

ggsave("Epi_curves_R0.jpeg", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)


### SIGMA EPI-CURVE PLOTTING ----
sig_varied <- all_data %>% filter(sigma != 0.5)

# preparing to plot
plot_list <- list() #storage of individual plots

efficacies <- unique(sig_varied$vaccine_efficacy)

# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- sig_varied %>% filter(vaccine_efficacy == v) %>%
    mutate(parm_cat = if_else(sigma > 0.5, "High", "Low"))
  
  # plot generation
  p <- ggplot()+
    geom_line(data = data_subset,
              aes(x = time, y = I, group=run_number, color = parm_cat),alpha = 0.1)+
    scale_color_manual(name = NULL, # Set to "Legend" if you want a title
                       values = c("High" = "darkred", "Low" = "steelblue")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = "Time (Days)",
         y = "Number of Cases")+
    coord_cartesian(ylim = c(0, max(sig_varied$I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
plot_annotation(
  title = 'Stochastic Infection Trajectories by Vaccine Efficacy (\u03C3)')

print(final_plot)

ggsave("Epi_curves_sig.jpeg", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)

### GAMMA EPI-CURVE PLOTTING ----
gam_varied <- all_data %>% filter(gamma != 0.2)

# preparing to plot
plot_list <- list() #storage of individual plots

efficacies <- unique(gam_varied$vaccine_efficacy)

# loop for each vaccine efficacy 
for (v in efficacies) {
  
  # filter for particular efficacy
  data_subset <- gam_varied %>% filter(vaccine_efficacy == v) %>%
    mutate(parm_cat = if_else(gamma > 0.2, "High", "Low"))
  
  # plot generation
  p <- ggplot()+
    geom_line(data = data_subset,
              aes(x = time, y = I, group=run_number, color = parm_cat),alpha = 0.1)+
    scale_color_manual(name = NULL, # Set to "Legend" if you want a title
                       values = c("High" = "darkred", "Low" = "steelblue")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    labs(title = paste0("Efficacy: ", v*100, "%"),
         x = "Time (Days)",
         y = "Number of Cases")+
    coord_cartesian(ylim = c(0, max(gam_varied$I)))+
    theme_minimal()
  
  # add to plot list
  plot_list[[as.character(v)]] <- p
  
}

final_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
plot_annotation(
  title = 'Stochastic Infection Trajectories by Vaccine Efficacy (\u03B3)')

print(final_plot)

ggsave("Epi_curves_gam.jpeg", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)
