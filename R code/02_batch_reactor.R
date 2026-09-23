# =============================================================================
# Batch Reactor Model — Paracetamol Synthesis
# =============================================================================
# Isothermal batch reactor with two reactions:
#   R1: A + B -> P + S  (desired, rate r1 = k1*CA*CB)
#   R2: A + B -> I + S  (undesired impurity, rate r2 = k2*CA*CB)
#
# Mole balances (V = const):
#   dCA/dt = -(r1 + r2)
#   dCB/dt = -(r1 + r2)   [same stoichiometry for both here]
#   dCP/dt =  r1
#   dCI/dt =  r2
# =============================================================================

source("R code/01_parameters.R")
library(deSolve)

batch_odes <- function(t, state, parms) {
  CA <- state["CA"]
  CB <- state["CB"]
  CP <- state["CP"]
  CI <- state["CI"]

  k1 <- k_arrhenius(parms$k1_ref, parms$Ea1, parms$T_batch)
  k2 <- k_arrhenius(parms$k2_ref, parms$Ea2, parms$T_batch)

  r1 <- k1 * CA * CB   # desired
  r2 <- k2 * CA * CB   # undesired

  list(c(
    dCA = -(r1 + r2),
    dCB = -(r1 + r2),
    dCP =  r1,
    dCI =  r2
  ))
}

run_batch <- function(params) {
  y0 <- c(CA = params$CA0, CB = params$CB0,
          CP = params$CP0, CI = params$CI0)

  times <- seq(0, params$t_batch_end, by = 0.5)

  out <- ode(y = y0, times = times, func = batch_odes,
             parms = params, method = "lsoda")
  df <- as.data.frame(out)

  # Derived quantities
  df$X_A      <- (params$CA0 - df$CA) / params$CA0       # conversion of A
  df$Y_P      <- df$CP / params$CA0                       # yield of paracetamol
  df$S_P      <- ifelse(df$CP + df$CI > 0,
                        df$CP / (df$CP + df$CI), NA)      # selectivity
  df$purity   <- df$CP / (df$CP + df$CI + 1e-9)          # mole fraction purity
  df$moles_P <- df$CP * params$V_batch_total     # total moles of paracetamol
  df$moles_I <- df$CI * params$V_batch_total     # impurity
  
  df$productivity <- ifelse(df$time > 0,
                            df$moles_P / df$time,
                            NA)
  df$reactor  <- "Batch"
  df
}

batch_results <- run_batch(params)

# Summary at end of batch
final <- tail(batch_results, 1)
cat("\n========== BATCH REACTOR RESULTS ==========\n")
cat(sprintf("  Reaction time        : %.1f min\n", max(batch_results$time)))
cat(sprintf("  Temperature          : %.0f K (%.0f°C)\n", params$T_batch, params$T_batch - 273))
cat(sprintf("  Final conversion X_A : %.3f (%.1f%%)\n", final$X_A, final$X_A*100))
cat(sprintf("  Paracetamol conc     : %.4f mol/L\n", final$CP))
cat(sprintf("  Impurity conc        : %.4f mol/L\n", final$CI))
cat(sprintf("  Selectivity S_P      : %.4f\n", final$S_P))
cat(sprintf("  Product purity       : %.2f%%\n", final$purity * 100))
cat(sprintf("  Yield Y_P            : %.4f\n", final$Y_P))
cat(sprintf("  Total paracetamol (mol) : %.2f\n", final$moles_P))
cat(sprintf("  Productivity (mol/min)  : %.2f\n", final$productivity))

# Time to reach target conversion
target_row <- batch_results[which(batch_results$X_A >= params$target_conv)[1], ]
if (!is.na(target_row$time)) {
  cat(sprintf("  Time to %.0f%% conv     : %.1f min\n",
              params$target_conv*100, target_row$time))
}

