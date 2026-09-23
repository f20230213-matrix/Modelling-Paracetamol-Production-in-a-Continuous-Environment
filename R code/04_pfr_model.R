# =============================================================================
# Plug Flow Reactor (PFR) Model — Paracetamol Synthesis
# =============================================================================
# PFR mole balances (plug flow, isothermal, V = axial coordinate):
#   dFA/dV = -r1 - r2          => dCA/d(V/Q) = -(r1+r2)
#   dFB/dV = -r1 - r2
#   dFP/dV =  r1
#   dFI/dV =  r2
#
# With tau = V/Q as the "space time" coordinate:
#   dCA/dtau = -(r1+r2)
#   dCB/dtau = -(r1+r2)
#   dCP/dtau =  r1
#   dCI/dtau =  r2
# (same ODE form as batch, but tau is spatial!)
# =============================================================================

source("R code/01_parameters.R")
library(deSolve)

pfr_odes <- function(tau, state, parms) {
  CA <- state["CA"]; CB <- state["CB"]
  CP <- state["CP"]; CI <- state["CI"]

  k1 <- k_arrhenius(parms$k1_ref, parms$Ea1, parms$T_pfr)
  k2 <- k_arrhenius(parms$k2_ref, parms$Ea2, parms$T_pfr)

  r1 <- k1 * CA * CB
  r2 <- k2 * CA * CB

  list(c(
    dCA = -(r1 + r2),
    dCB = -(r1 + r2),
    dCP =  r1,
    dCI =  r2
  ))
}

run_pfr <- function(params) {
  y0 <- c(CA = params$CA0, CB = params$CB0,
          CP = params$CP0, CI = params$CI0)

  tau_max <- params$V_pfr_total / params$Q_pfr   # total space time [min]
  taus <- seq(0, tau_max, length.out = 500)

  out <- ode(y = y0, times = taus, func = pfr_odes,
             parms = params, method = "lsoda")
  df <- as.data.frame(out)
  names(df)[1] <- "tau"

  df$V_axial <- df$tau * params$Q_pfr   # axial volume coordinate [L]
  df$X_A     <- (params$CA0 - df$CA) / params$CA0
  df$Y_P     <- df$CP / params$CA0
  df$S_P     <- df$CP / (df$CP + df$CI + 1e-9)
  df$purity  <- df$CP / (df$CP + df$CI + 1e-9)
  df$production_rate <- params$Q_pfr * df$CP
  df$reactor <- "PFR"
  df
}

# Volume sweep — how much PFR volume for target conversion?
run_pfr_volume_sweep <- function(params) {
  Q     <- params$Q_pfr
  tau_max <- 40   # min maximum space time to sweep
  taus  <- seq(0.1, tau_max, by = 0.5)

  y0 <- c(CA = params$CA0, CB = params$CB0,
          CP = params$CP0, CI = params$CI0)
  out <- ode(y = y0, times = c(0, taus), func = pfr_odes,
             parms = params, method = "lsoda")
  df <- as.data.frame(out)
  names(df)[1] <- "tau"
  df$V_reactor <- df$tau * Q
  df$X_A   <- (params$CA0 - df$CA) / params$CA0
  df$Y_P   <- df$CP / params$CA0
  df$S_P   <- df$CP / (df$CP + df$CI + 1e-9)
  df$purity <- df$CP / (df$CP + df$CI + 1e-9)
  df$production_rate <- Q * df$CP
  df
}

pfr_results <- run_pfr(params)
pfr_sweep   <- run_pfr_volume_sweep(params)

final_pfr <- tail(pfr_results, 1)
tau_total  <- params$V_pfr_total / params$Q_pfr

cat("\n========== PFR RESULTS ==========\n")
cat(sprintf("  Reactor volume       : %.0f L\n",    params$V_pfr_total))
cat(sprintf("  Volumetric flow rate : %.1f L/min\n", params$Q_pfr))
cat(sprintf("  Space time tau       : %.2f min\n",  tau_total))
cat(sprintf("  Temperature          : %.0f K (%.0f°C)\n", params$T_pfr, params$T_pfr - 273))
cat(sprintf("  Conversion X_A       : %.3f (%.1f%%)\n", final_pfr$X_A, final_pfr$X_A*100))
cat(sprintf("  Paracetamol conc     : %.4f mol/L\n", final_pfr$CP))
cat(sprintf("  Impurity conc        : %.4f mol/L\n", final_pfr$CI))
cat(sprintf("  Selectivity S_P      : %.4f\n",       final_pfr$S_P))
cat(sprintf("  Product purity       : %.2f%%\n",     final_pfr$purity * 100))
cat(sprintf("  Yield Y_P            : %.4f\n",       final_pfr$Y_P))
cat(sprintf("  Productivity (mol/min)  : %.2f\n", final_pfr$production_rate))

# Space time / volume needed for 90% conversion
target_row <- pfr_sweep[which(pfr_sweep$X_A >= params$target_conv)[1], ]
if (!is.na(target_row$tau)) {
  cat(sprintf("  Tau for 90%% conv     : %.2f min\n", target_row$tau))
  cat(sprintf("  Volume for 90%% conv  : %.1f L\n",   target_row$V_reactor))
}
  
