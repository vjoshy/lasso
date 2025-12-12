# ==============================================================================
# Description: Simulation results analysis script - plotting and summary
# ==============================================================================

library(dplyr)

final_simulation_results <- readRDS("data/final_simulation_results.rds") %>%
  mutate(Method_Name = ifelse(Method_Name == "Ground Truth", "glmnet", Method_Name))


# final_simulation_results %>%
#   group_by(p, n, Method_Type) %>%
#   summarize(
#     precision = mean(Precision), 
#     recall = mean(Recall), 
#     f1 = mean(F1_Score),
#     .groups = "drop"
#   ) %>%
#   View()


library(tidyr)
library(ggplot2)

# Summarize results by Scenario and Method
simulation_summary <- final_simulation_results %>%
  group_by(n, p, sparsity, corr_type, snr, Method_Name) %>%
  summarize(
    mean_Precision = mean(Precision, na.rm = TRUE),
    sd_Precision   = sd(Precision, na.rm = TRUE),
    mean_Recall    = mean(Recall, na.rm = TRUE),
    sd_Recall      = sd(Recall, na.rm = TRUE),
    mean_F1_Score  = mean(F1_Score, na.rm = TRUE),
    sd_F1_Score    = sd(F1_Score, na.rm = TRUE),
    mean_Runtime   = mean(Runtime, na.rm = TRUE),
    sd_Runtime     = sd(Runtime, na.rm = TRUE),
    mean_MSE       = mean(MSE, na.rm = TRUE),
    
    .groups = "drop"
  )

print(head(simulation_summary))


# Summarize by Dimension (p) and Method
table_data <- final_simulation_results %>%
  group_by(p, Method_Name) %>%
  summarize(
    mse = mean(MSE),
    MSE_SD = sd(MSE, na.rm = TRUE),
    run = mean(Runtime),
    sd_rtun = sd(Runtime, na.rm = TRUE),
    F1 = mean(F1_Score),
    sd_f1 = sd(F1_Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p, Method_Name)

plot_data_stats <- simulation_summary %>%
  pivot_longer(
    cols = c(mean_Precision, mean_Recall, mean_F1_Score),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  mutate(
    Metric = gsub("mean_", "", Metric),
    Metric = gsub("_", " ", Metric),
    p_label = as.factor(p),
    n_label = paste0("n = ", n)
  )

# Define nice labels 
metric_labels <- c(
  "Precision" = "Precision",
  "Recall"    = "Recall",
  "F1 Score"  = "F1 Score"
)


ggplot(plot_data_stats, aes(x = p_label, y = Score, fill = Method_Name)) +

  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(Metric ~ n_label, labeller = labeller(Metric = as_labeller(metric_labels))) +

  scale_fill_brewer(palette = "Set1", name = "Method") +
  coord_cartesian(ylim = c(0, 1.05)) + 
  theme_bw() +
  labs(
    title = "Statistical Performance Comparison",
    subtitle = "Comparison of Precision, Recall, and F1 across dimensions (p) and sample sizes (n)",
    y = "Mean Score",
    x = "Number of Covariates (p)"
  ) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


ggplot(simulation_summary, aes(x = p, y = mean_Runtime, color = Method_Name)) +
  
  geom_line(linewidth = 1) + 
  geom_point(size = 3) +
  
  facet_grid(n ~ sparsity, labeller = label_both) +
  
  scale_color_brewer(palette = "Set1", name = "Method") +
  scale_x_continuous(breaks = c(10, 100, 200)) + 
  scale_y_log10() + 
  theme_bw() +
  labs(
    title = "Computational Efficiency: Runtime vs. Dimension",
    subtitle = "Impact of Screening Rules across sparsity levels (Log Scale Y-axis)",
    y = "Mean Runtime (seconds)",
    x = "Dimension (p)"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )


# 1. Prepare Data
plot_data_methods <- simulation_summary %>%
  mutate(
    snr_label = paste("SNR =", snr),
    p_label = paste("p =", p)
  )

ggplot(plot_data_methods, aes(x = as.factor(n), y = mean_F1_Score, fill = as.factor(snr))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  
  facet_grid(Method_Name ~ p_label) +
  
  scale_fill_brewer(palette = "Paired", name = "Signal-to-Noise Ratio") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  labs(
    title = "Impact of SNR on Variable Selection by Method",
    subtitle = "Comparison across Sample Sizes (n) and Dimensions (p)",
    y = "Mean F1-Score",
    x = "Sample size (n)"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggplot(plot_data_methods, aes(x = corr_type, y = mean_F1_Score, fill = corr_type)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
    facet_grid(Method_Name ~ p, labeller = labeller(p = label_both)) +

  scale_fill_manual(values = c("independent" = "#1f77b4", "block_diagonal" = "#d62728"),
                    labels = c("Block Diagonal", "Independent"),
                    name = "Correlation Structure") +
  scale_x_discrete(labels = c("Block", "Indep")) +
  theme_bw() +
  labs(
    title = "Impact of Feature Correlation by Method",
    y = "Mean F1-Score",
    x = "Correlation Structure"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )
