library(dplyr)
library(ggplot2)
library(cowplot)
library(tornado)

# Load econ model results for Single Value Cost Estimate ----
econ_results <- readRDS("Baseline_Unique_Outcomes_Single_Est.rds")

# Summary Stats for Single Value Cost Estimate ----
Summary_stats <- econ_results %>% 
  group_by(VE, vaccine_cost_per_pig) %>% 
  summarise(
    
    # sample size
    N = n(),
    
    #NPV
    NPV_mean = mean(NPV),
    NPV_sd = sd(NPV),
    NPV_median = median(NPV),
    
    #BCR
    BCR_mean = mean(BCR),
    BCR_sd = sd(BCR),
    BCR_median = median(BCR),
    
    #ROI
    ROI_mean = mean(ROI),
    ROI_sd = sd(ROI),
    ROI_median = median(ROI),
    
    .groups = "drop"
  ) %>%
  mutate(U_CI_NPV = NPV_mean + 1.96 * NPV_sd/sqrt(N),
         L_CI_NPV = NPV_mean - 1.96 * NPV_sd/sqrt(N),
         U_CI_BCR = BCR_mean + 1.96 * BCR_sd/sqrt(N),
         L_CI_BCR = BCR_mean - 1.96 * BCR_sd/sqrt(N),
         U_CI_ROI = ROI_mean + 1.96 * ROI_sd/sqrt(N),
         L_CI_ROI = ROI_mean - 1.96 * ROI_sd/sqrt(N)) %>%
  mutate(PP_NPV_mean)
  relocate( VE, vaccine_cost_per_pig,
            NPV_mean, NPV_sd, NPV_median, L_CI_NPV, U_CI_NPV, 
            BCR_mean, BCR_sd, BCR_median, L_CI_BCR, U_CI_BCR,
            ROI_mean, ROI_sd, ROI_median, L_CI_ROI, U_CI_ROI)

write.csv(Summary_stats, "Baseline_Single_est_sum_stat.csv", row.names = F)

# Load econ MC model for varied 
MC_econ_results <- readRDS("Baseline_Unique_Outcomes_MC.rds") 

# Prop of BCR > 1 
VE <- unique(MC_econ_results$VE)
vac_cost <- unique(MC_econ_results$vaccine_cost_per_pig)

# setting up the dataframe
MC_BCR <- data.frame()

MC_BCR <- MC_econ_results %>%
  group_by(VE, vaccine_cost_per_pig) %>%
  summarise(
    count_BCR_gt1 = sum(BCR >= 1, na.rm = T),
    total_observations = n(),
    mean_bcr = mean(BCR),
    BCR_sd = sd(BCR)
  ) %>% 
  mutate(percentage = count_BCR_gt1/total_observations * 100) %>%
  mutate(U_CI_BCR = mean_bcr + 1.96 * BCR_sd/sqrt(total_observations),
         L_CI_BCR = mean_bcr - 1.96 * BCR_sd/sqrt(total_observations),
         )
write.csv(MC_BCR, "BCR_MC_stat.csv", row.names = F)

# Prop of ROI > 10%
MC_ROI <- data.frame()
MC_ROI <- MC_econ_results %>%
  filter(BCR > 1) %>%
  group_by(VE, vaccine_cost_per_pig) %>%
  summarise(
    count_gt5 = sum(ROI >= 10, na.rm = T ),
    total_obs = n()
  ) %>%
  mutate(percent = count_gt5/total_obs *100)

average_ROI <- mean(MC_ROI$percent)

ROI_stats <- data.frame()
ROI_stats <- MC_econ_results %>%
  group_by(VE, vaccine_cost_per_pig) %>% 
  summarise(
    
    # sample size
    N = n(),
    
    #ROI
    ROI_mean = mean(ROI),
    ROI_sd = sd(ROI),
    ROI_median = median(ROI),
    count_gt10 = sum(ROI >= 10, na.rm = T ),
    
    .groups = "drop"
  ) %>%
  mutate(percent = count_gt10/N *100)
  

# Make sure VE and vaccine cost are treated as factors for plotting
MC_BCR <- MC_BCR%>%
  mutate(
    VE = factor(VE),
    vaccine_cost_per_pig = factor(vaccine_cost_per_pig)
  )


# Create the heat map - PERCENTAGE GREATER THAN 1
p1_heat<-ggplot(MC_BCR, aes(x = VE, y = vaccine_cost_per_pig, fill = percentage)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "navy", 
                       midpoint = 50, limits = c(0, 100),
                       name = "% BCR >= 1") +
  labs(x = "Vaccine Efficacy",y = "Vaccine Cost per Pig ($)",) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

print(p1_heat)

#ggsave("HM_BCR_MC_gt1.jpeg",device = "jpeg", plot = p1_heat, width = 6, height = 4, units = "in", dpi = 600)


# Create the heat map - BCR Value
p2_heat<-ggplot(MC_BCR, aes(x = VE, y = vaccine_cost_per_pig, fill = mean_bcr)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "tomato", mid = "white", high = "navy", 
                       midpoint = 1, limits = c(0, 4.5),
                       name = "Mean BCR") +
  labs(x = "Vaccine Efficacy",y = "Vaccine Cost per Pig ($)",) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

print(p2_heat)

#ggsave("HM_BCR_MC_val.jpeg",device = "jpeg", plot = p2_heat, width = 6, height = 4, units = "in", dpi = 600)

#Merging BCR Plots together
cp <- plot_grid(p1_heat, p2_heat, labels = c('A', 'B'), ncol = 1, align = 'h', axis = "lr")

print(cp)

ggsave("BCR_Heatmap.jpeg", plot = cp, width = 6, height = 8, units = "in", dpi = 600)  
  
  
  
#### Tornado Plot ---- 
# Combining epi uncertainty and economic uncertainty
epi_econ <- readRDS("Uncertainty_Analysis_MC_SP.rds")

all_econ <- bind_rows(MC_econ_results, epi_econ)

# Getting tornado data from model
lm1 <- lm(BCR ~ VE+vaccine_cost_per_pig+penalty_kg+mortality+raising_cost+profit_loss+carcass_disposal+feed_cost_kg, 
          data = all_econ) #+
            #R0_Val+Sigma_Val+Gamma_Val, 
         # data = all_econ)

torn <- tornado::tornado(lm1, type = "ranges", alpha = 0.05)

# Basic plot
plot(torn, xlabel = "BCR", 
     geom_bar_control = list(width =0.4))


# Convert the tornado object to a data frame for ggplot2
# The 'torn' object has a 'data' component that's suitable for plotting
plot_data <- as.data.frame(torn$data)

plot_data <- plot_data %>%
  mutate(plotdat.variable = recode(plotdat.variable,
                                   "VE" = "Vaccine Efficacy",
                                   "vaccine_cost_per_pig" = "Vaccination Cost",
                                   "penalty_kg" = "Feed Penalty",
                                   "mortality" = "Mortality Rate",
                                   "raising_cost" = "Raising Pig Cost",
                                   "profit_loss" = "Lost Profit from Death",
                                   "carcass_disposal" = "Carcass Disposal",
                                   "feed_cost_kg" = "Feed Cost per kg" )) #,
                                   # "R0_Val" = "R0",
                                   # "Gamma_Val" = "\u03B3",
                                   # "Sigma_Val" = "\u03C3"))

# Create the ggplot USE THIS ONE RN
p_tornado <- ggplot(plot_data, aes(x = plotdat.value, y = plotdat.variable)) +
  geom_bar(stat = "identity", aes(fill = plotdat.Level), width = 0.4) + # Use geom_bar for bar plot
  scale_x_continuous(name = "Change in Benefit-Cost Ratio", labels = scales::number_format(accuracy = 0.01)) + # Format x-axis to 0.01 precision
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  labs(fill="Level") + 
  theme_minimal() + # Use a clean theme
  theme(axis.title.y = element_blank()) # Remove y-axis title as variables are clear

print(p_tornado)

ggsave("Tornado.tiff", plot = p_tornado, width = 6, height = 4, units = "in", dpi = 600)

  
  
  
  