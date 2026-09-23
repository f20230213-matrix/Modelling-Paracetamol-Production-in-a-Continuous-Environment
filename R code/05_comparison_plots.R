# =============================================================================
# Comparison & Visualisation — Paracetamol Reactor Study
# =============================================================================

source("R code/02_batch_reactor.R")
source("R code/03_cstr_model.R")
source("R code/04_pfr_model.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

theme_paracetamol <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
      legend.position = "bottom",
      legend.title  = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      plot.margin = margin(8, 10, 8, 10)
    )
}

pal <- c("Batch" = "#2E86AB", "CSTR" = "#E84855", "PFR" = "#3BB273")

cat("\n--- Generating plots ---\n")

# ===========================================================================
# PLOT 1: Concentration profiles (all three reactors)
# ===========================================================================

# Batch: time axis
batch_long <- batch_results %>%
  select(time, CA, CB, CP, CI, reactor) %>%
  pivot_longer(cols = c(CA, CB, CP, CI), names_to = "species",
               values_to = "concentration") %>%
  rename(x_axis = time) %>%
  mutate(x_label = "Time (min)")

# PFR: tau axis
pfr_long <- pfr_results %>%
  select(tau, CA, CB, CP, CI, reactor) %>%
  pivot_longer(cols = c(CA, CB, CP, CI), names_to = "species",
               values_to = "concentration") %>%
  rename(x_axis = tau) %>%
  mutate(x_label = "Space Time / min")

species_labels <- c(CA = "4-Aminophenol (A)",
                    CB = "Acetic Anhydride (B)",
                    CP = "Paracetamol (P)",
                    CI = "Impurity (I)")

p1_batch <- batch_long %>%
  mutate(species = factor(species, levels = c("CA","CB","CP","CI"))) %>%
  ggplot(aes(x = x_axis, y = concentration, color = species, linetype = species)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(CA="#264653", CB="#2A9D8F", CP="#E76F51", CI="#F4A261"),
                     labels = species_labels) +
  scale_linetype_manual(values = c(CA="solid", CB="dashed", CP="solid", CI="dotted"),
                        labels = species_labels) +
  labs(title = "Batch Reactor", subtitle = "T = 323 K (50°C)",
       x = "Time (min)", y = "Concentration (mol/L)") +
  ylim(0, 1.3) + theme_paracetamol()

# CSTR: horizontal lines at steady-state + dynamic approach
p1_cstr <- cstr_dynamic %>%
  select(time, CA, CB, CP, CI) %>%
  pivot_longer(cols = c(CA, CB, CP, CI), names_to = "species",
               values_to = "concentration") %>%
  mutate(species = factor(species, levels = c("CA","CB","CP","CI"))) %>%
  ggplot(aes(x = time, y = concentration, color = species, linetype = species)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(CA="#264653", CB="#2A9D8F", CP="#E76F51", CI="#F4A261"),
                     labels = species_labels) +
  scale_linetype_manual(values = c(CA="solid", CB="dashed", CP="solid", CI="dotted"),
                        labels = species_labels) +
  labs(title = "CSTR (Dynamic Startup)", subtitle = "T = 333 K (60°C), τ = 30 min",
       x = "Time (min)", y = "Concentration (mol/L)") +
  ylim(0, 1.3) + theme_paracetamol()

p1_pfr <- pfr_long %>%
  mutate(species = factor(species, levels = c("CA","CB","CP","CI"))) %>%
  ggplot(aes(x = x_axis, y = concentration, color = species, linetype = species)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(CA="#264653", CB="#2A9D8F", CP="#E76F51", CI="#F4A261"),
                     labels = species_labels) +
  scale_linetype_manual(values = c(CA="solid", CB="dashed", CP="solid", CI="dotted"),
                        labels = species_labels) +
  labs(title = "PFR", subtitle = "T = 328 K (55°C), V = 100 L",
       x = "Space Time (min)", y = "Concentration (mol/L)") +
  ylim(0, 1.3) + theme_paracetamol()

g1 <- arrangeGrob(p1_batch, p1_cstr, p1_pfr, ncol = 3,
                  top = grid::textGrob("Concentration Profiles — Paracetamol Synthesis",
                                       gp = grid::gpar(fontsize = 14, fontface = "bold")))
ggsave("Plots/01_concentration_profiles.png", g1, width = 15, height = 5, dpi = 180, type = "cairo")
cat("  Saved: plots/01_concentration_profiles.png\n")

# ===========================================================================
# PLOT 2: Conversion & Yield comparison
# ===========================================================================

# Unified time/tau axis for comparison
batch_kpi <- batch_results %>%
  select(time, X_A, Y_P, S_P, purity, productivity, reactor) %>%
  rename(x = time, production_rate = productivity)

pfr_kpi <- pfr_results %>%
  select(tau, X_A, Y_P, S_P, purity, production_rate, reactor) %>%
  rename(x = tau)

# For CSTR: show evolution during startup
cstr_kpi <- cstr_dynamic %>%
  select(time, X_A, Y_P, S_P, purity, production_rate, reactor) %>%
  rename(x = time)

all_kpi <- bind_rows(batch_kpi, pfr_kpi, cstr_kpi) %>%
  mutate(reactor = factor(reactor, levels = c("Batch","CSTR","PFR")))

p2a <- ggplot(all_kpi, aes(x = x, y = X_A, color = reactor)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  annotate("text", x = 5, y = 0.92, label = "90% target", size = 3, color = "grey40") +
  scale_color_manual(values = pal) +
  labs(title = "Conversion of 4-Aminophenol",
       x = "Time / Space Time (min)", y = "Conversion X_A") +
  theme_paracetamol()

p2b <- ggplot(all_kpi, aes(x = x, y = Y_P, color = reactor)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pal) +
  labs(title = "Yield of Paracetamol",
       x = "Time / Space Time (min)", y = "Yield Y_P") +
  theme_paracetamol()

p2c <- ggplot(all_kpi, aes(x = x, y = S_P, color = reactor)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pal) +
  labs(title = "Selectivity to Paracetamol",
       x = "Time / Space Time (min)", y = "Selectivity S_P") +
  theme_paracetamol()

p2d <- ggplot(all_kpi, aes(x = x, y = purity * 100, color = reactor)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pal) +
  labs(title = "Product Purity",
       x = "Time / Space Time (min)", y = "Purity (mol%)") +
  theme_paracetamol()

p_prod <- ggplot(all_kpi, aes(x = x, y = production_rate, color = reactor)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pal) +
  labs(title = "Production Rate Comparison",
       x = "Time / Space Time (min)", y = "Production Rate (mol/min)") +
  theme_paracetamol()

g2 <- arrangeGrob(p2a, p2b, p2c, p2d, p_prod, ncol = 3,
                  top = grid::textGrob("Performance Metrics — Reactor Comparison",
                                       gp = grid::gpar(fontsize = 14, fontface = "bold")))
ggsave("Plots/02_performance_metrics.png", g2, width = 12, height = 9, dpi = 180, type = "cairo")
cat("  Saved: plots/02_performance_metrics.png\n")

# ===========================================================================
# PLOT 3: Residence time / space time sensitivity
# ===========================================================================

# CSTR tau sweep
p3a <- ggplot(cstr_sweep, aes(x = tau)) +
  geom_line(aes(y = X_A, color = "Conversion"), linewidth = 1.2) +
  geom_line(aes(y = Y_P, color = "Yield"), linewidth = 1.2) +
  geom_line(aes(y = purity, color = "Purity"), linewidth = 1.2) +
  geom_vline(xintercept = params$tau_cstr, linetype = "dashed", color = "grey40") +
  annotate("text", x = params$tau_cstr + 2, y = 0.2,
           label = sprintf("τ = %d min", params$tau_cstr), size = 3, color = "grey30") +
  scale_color_manual(values = c(Conversion="#2E86AB", Yield="#3BB273", Purity="#E84855")) +
  labs(title = "CSTR: Effect of Residence Time",
       x = "Residence Time τ (min)", y = "Fraction") +
  theme_paracetamol()

# PFR volume sweep
p3b <- ggplot(pfr_sweep, aes(x = V_reactor)) +
  geom_line(aes(y = X_A, color = "Conversion"), linewidth = 1.2) +
  geom_line(aes(y = Y_P, color = "Yield"), linewidth = 1.2) +
  geom_line(aes(y = purity, color = "Purity"), linewidth = 1.2) +
  geom_vline(xintercept = params$V_pfr_total, linetype = "dashed", color = "grey40") +
  annotate("text", x = params$V_pfr_total + 3, y = 0.2,
           label = sprintf("V = %d L", params$V_pfr_total), size = 3, color = "grey30") +
  scale_color_manual(values = c(Conversion="#2E86AB", Yield="#3BB273", Purity="#E84855")) +
  labs(title = "PFR: Effect of Reactor Volume",
       x = "Reactor Volume (L)", y = "Fraction") +
  theme_paracetamol()

g3 <- arrangeGrob(p3a, p3b, ncol = 2,
                  top = grid::textGrob("Sensitivity Analysis — CSTR τ and PFR Volume",
                                       gp = grid::gpar(fontsize = 14, fontface = "bold")))
ggsave("Plots/03_sensitivity_analysis.png", g3, width = 12, height = 5, dpi = 180, type = "cairo")
cat("  Saved: plots/03_sensitivity_analysis.png\n")

# ===========================================================================
# PLOT 4: Temperature effect on selectivity (k1/k2 ratio)
# ===========================================================================
temp_range <- seq(290, 370, by = 2)
temp_df <- data.frame(
  T = temp_range,
  k1 = sapply(temp_range, function(T) k_arrhenius(params$k1_ref, params$Ea1, T)),
  k2 = sapply(temp_range, function(T) k_arrhenius(params$k2_ref, params$Ea2, T))
) %>% mutate(
  ratio = k1 / k2,
  S_intrinsic = k1 / (k1 + k2)
)

p4a <- ggplot(temp_df, aes(x = T - 273)) +
  geom_line(aes(y = k1, color = "k1 (desired)"), linewidth = 1.2) +
  geom_line(aes(y = k2, color = "k2 (undesired)"), linewidth = 1.2) +
  geom_vline(xintercept = c(params$T_batch, params$T_cstr, params$T_pfr) - 273,
             linetype = "dotted", color = c("#2E86AB","#E84855","#3BB273"), linewidth = 0.8) +
  scale_color_manual(values = c("k1 (desired)" = "#3BB273", "k2 (undesired)" = "#E84855")) +
  labs(title = "Rate Constants vs Temperature",
       x = "Temperature (°C)", y = "k (L/mol·min)") +
  theme_paracetamol()

p4b <- ggplot(temp_df, aes(x = T - 273, y = S_intrinsic * 100)) +
  geom_line(linewidth = 1.3, color = "#2E86AB") +
  geom_vline(xintercept = c(params$T_batch, params$T_cstr, params$T_pfr) - 273,
             linetype = "dotted", color = c("#2E86AB","#E84855","#3BB273"), linewidth = 0.8) +
  annotate("text", x = params$T_batch-273+1.5, y = 93.5,
           label = "Batch\n50°C", size = 2.8, color = "#2E86AB") +
  annotate("text", x = params$T_cstr-273+1.5, y = 93.5,
           label = "CSTR\n60°C", size = 2.8, color = "#E84855") +
  annotate("text", x = params$T_pfr-273+1.5, y = 93.5,
           label = "PFR\n55°C", size = 2.8, color = "#3BB273") +
  labs(title = "Intrinsic Selectivity vs Temperature",
       x = "Temperature (°C)", y = "Selectivity S_P (%)") +
  theme_paracetamol()

g4 <- arrangeGrob(p4a, p4b, ncol = 2,
                  top = grid::textGrob("Temperature Effects on Kinetics and Selectivity",
                                       gp = grid::gpar(fontsize = 14, fontface = "bold")))
ggsave("Plots/04_temperature_effects.png", g4, width = 12, height = 5, dpi = 180, type = "cairo")
cat("  Saved: plots/04_temperature_effects.png\n")

# ===========================================================================
# PLOT 5: Summary bar chart comparison
# ===========================================================================
final_batch <- tail(batch_results, 1)
final_cstr  <- cstr_ss
final_pfr   <- tail(pfr_results, 1)

summary_df <- data.frame(
  Reactor   = c("Batch", "CSTR", "PFR"),
  Conversion = c(final_batch$X_A, final_cstr$X_A, final_pfr$X_A),
  Yield      = c(final_batch$Y_P, final_cstr$Y_P, final_pfr$Y_P),
  Selectivity= c(final_batch$S_P, final_cstr$S_P, final_pfr$S_P),
  Purity     = c(final_batch$purity, final_cstr$purity, final_pfr$purity)
) %>%
  pivot_longer(-Reactor, names_to = "Metric", values_to = "Value") %>%
  mutate(Reactor = factor(Reactor, levels = c("Batch","CSTR","PFR")),
         Metric  = factor(Metric, levels = c("Conversion","Yield","Selectivity","Purity")))

p5 <- ggplot(summary_df, aes(x = Metric, y = Value * 100, fill = Reactor)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", Value * 100)),
            position = position_dodge(0.7), vjust = -0.4, size = 3) +
  scale_fill_manual(values = pal) +
  labs(title = "Reactor Performance Comparison — Paracetamol Synthesis",
       subtitle = "Batch (60 min) vs CSTR (τ=30 min) vs PFR (V=100 L, Q=3 L/min)",
       x = NULL, y = "Value (%)") +
  theme_paracetamol() +
  theme(legend.position = "right")

ggsave("Plots/05_summary_comparison.png", p5, width = 10, height = 6, dpi = 180, type = "cairo")
cat("  Saved: plots/05_summary_comparison.png\n")

cat("\n--- All plots generated successfully. ---\n")

