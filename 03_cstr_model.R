# =============================================================================
# Continuous CSTR Model — Paracetamol Synthesis
# =============================================================================
# Steady-state CSTR balance (perfect mixing, V = const, isothermal):
#   0 = CA0 - CA - tau*(r1+r2)
#   0 = CB0 - CB - tau*(r1+r2)
#   0 = -CP  + tau*r1
#   0 = -CI  + tau*r2
# where tau = V/Q (residence time)
#
# Also solve startup transient (dynamic CSTR) as ODE:
#   V*dCi/dt = Q*(Ci_feed - Ci) + V*ri
# =============================================================================

source("C:/Users/inbox/Downloads/paracetamol_reactor_study/paracetamol_project/R/01_parameters.R")
library(deSolve)

# ---------- Steady-state (algebraic solve via iteration) ----------
cstr_steady_state <- function(params, tau = params$tau_cstr) {
  k1 <- k_arrhenius(params$k1_ref, params$Ea1, params$T_cstr)
  k2 <- k_arrhenius(params$k2_ref, params$Ea2, params$T_cstr)

  CA0 <- params$CA0; CB0 <- params$CB0

  # For A + B -> products: r = k*CA*CB
  # CSTR: CA = CA0 - tau*(k1+k2)*CA*CB
  #        CB = CB0 - tau*(k1+k2)*CA*CB
  # So CA0 - CA = CB0 - CB  => CB = CA + (CB0 - CA0)
  # Substituting: CA = CA0 - tau*(k1+k2)*CA*(CA + CB0-CA0)
  # Solve as quadratic in CA

  delta <- CB0 - CA0
  k_tot <- (k1 + k2) * tau

  # CA * (1 + k_tot*(CA + delta)) = CA0
  # k_tot*CA^2 + (1 + k_tot*delta)*CA - CA0 = 0
  a_q <- k_tot
  b_q <- 1 + k_tot * delta
  c_q <- -CA0
  disc <- b_q^2 - 4*a_q*c_q
  CA_ss <- (-b_q + sqrt(disc)) / (2 * a_q)
  CB_ss <- CA_ss + delta

  r1_ss <- k1 * CA_ss * CB_ss
  r2_ss <- k2 * CA_ss * CB_ss

  CP_ss <- tau * r1_ss
  CI_ss <- tau * r2_ss
  # Production rate (mol/min)

  list(CA = CA_ss, CB = CB_ss, CP = CP_ss, CI = CI_ss,
       X_A = (CA0 - CA_ss)/CA0,
       S_P = CP_ss / (CP_ss + CI_ss),
       purity = CP_ss / (CP_ss + CI_ss),
       Y_P = CP_ss / CA0,
       production_rate = params$Q_pfr * CP_ss,
       tau = tau, k1 = k1, k2 = k2)
}

# ---------- Dynamic startup (ODE) ----------
cstr_dynamic_odes <- function(t, state, parms) {
  CA <- state["CA"]; CB <- state["CB"]
  CP <- state["CP"]; CI <- state["CI"]

  k1 <- k_arrhenius(parms$k1_ref, parms$Ea1, parms$T_cstr)
  k2 <- k_arrhenius(parms$k2_ref, parms$Ea2, parms$T_cstr)

  r1 <- k1 * CA * CB
  r2 <- k2 * CA * CB

  tau <- parms$tau_cstr

  list(c(
    dCA = (parms$CA0 - CA)/tau - (r1 + r2),
    dCB = (parms$CB0 - CB)/tau - (r1 + r2),
    dCP = (0         - CP)/tau + r1,
    dCI = (0         - CI)/tau + r2
  ))
}

run_cstr_dynamic <- function(params) {
  # Start from empty reactor
  y0 <- c(CA = params$CA0, CB = params$CB0, CP = 0, CI = 0)
  times <- seq(0, 5 * params$tau_cstr, by = 0.5)

  out <- ode(y = y0, times = times, func = cstr_dynamic_odes,
             parms = params, method = "lsoda")
  df <- as.data.frame(out)
  df$X_A    <- (params$CA0 - df$CA) / params$CA0
  df$Y_P    <- df$CP / params$CA0
  df$S_P    <- df$CP / (df$CP + df$CI + 1e-9)
  df$purity <- df$CP / (df$CP + df$CI + 1e-9)
  df$production_rate <- params$Q_pfr * df$CP
  df$reactor <- "CSTR"
  df
}

# ---------- Tau sensitivity analysis ----------
run_cstr_tau_sweep <- function(params) {
  tau_vals <- seq(5, 120, by = 5)
  results <- lapply(tau_vals, function(tau) {
    ss <- cstr_steady_state(params, tau)
    data.frame(tau = tau, X_A = ss$X_A, Y_P = ss$Y_P,
               S_P = ss$S_P, purity = ss$purity,
               CP = ss$CP, CI = ss$CI,
               production_rate = params$Q_pfr * ss$CP)
  })
  do.call(rbind, results)
}

# Run everything
cstr_ss      <- cstr_steady_state(params)
cstr_dynamic <- run_cstr_dynamic(params)
cstr_sweep   <- run_cstr_tau_sweep(params)

cat("\n========== CSTR STEADY-STATE RESULTS ==========\n")
cat(sprintf("  Residence time tau   : %.1f min\n",   params$tau_cstr))
cat(sprintf("  Temperature          : %.0f K (%.0f°C)\n", params$T_cstr, params$T_cstr - 273))
cat(sprintf("  Conversion X_A       : %.3f (%.1f%%)\n", cstr_ss$X_A, cstr_ss$X_A*100))
cat(sprintf("  Paracetamol conc     : %.4f mol/L\n", cstr_ss$CP))
cat(sprintf("  Impurity conc        : %.4f mol/L\n", cstr_ss$CI))
cat(sprintf("  Selectivity S_P      : %.4f\n",       cstr_ss$S_P))
cat(sprintf("  Product purity       : %.2f%%\n",     cstr_ss$purity * 100))
cat(sprintf("  Yield Y_P            : %.4f\n",       cstr_ss$Y_P))
cat(sprintf("  Productivity (mol/min)  : %.2f\n", cstr_ss$production_rate))

# Tau required for 90% conversion
tau_90 <- cstr_sweep$tau[which(cstr_sweep$X_A >= params$target_conv)[1]]
cat(sprintf("  Tau for 90%% conv     : %.0f min\n",  tau_90))

