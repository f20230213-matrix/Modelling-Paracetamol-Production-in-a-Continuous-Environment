# =============================================================================
# Paracetamol Reactor Study — Master Run Script
# Generates all results, tables, and the final comparison report
# =============================================================================


cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   PARACETAMOL SYNTHESIS — REACTOR MODELLING STUDY               ║\n")
cat("║   Batch vs Continuous CSTR vs PFR                               ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# Run all models and generate plots
source("R code/05_comparison_plots.R")

# ===========================================================================
# FINAL SUMMARY TABLE
# ===========================================================================
final_batch <- tail(batch_results, 1)
final_pfr   <- tail(pfr_results, 1)

tau_pfr <- params$V_pfr_total / params$Q_pfr

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║                FINAL PERFORMANCE SUMMARY TABLE                         ║\n")
cat("╠══════════════════════════════════════════════╦═════════╦══════════╦══════╣\n")
cat("║ Metric                                       ║  Batch  ║   CSTR   ║  PFR ║\n")
cat("╠══════════════════════════════════════════════╬═════════╬══════════╬══════╣\n")
cat(sprintf("║ Temperature (°C)                             ║   %3.0f   ║    %3.0f   ║  %3.0f ║\n",
    params$T_batch-273, params$T_cstr-273, params$T_pfr-273))
cat(sprintf("║ Reaction time / space time (min)             ║   %3.0f   ║    %3.0f   ║ %4.1f ║\n",
    params$t_batch_end, params$tau_cstr, tau_pfr))
cat(sprintf("║ Conversion X_A (%%)                           ║  %5.1f  ║   %5.1f  ║ %4.1f ║\n",
    final_batch$X_A*100, cstr_ss$X_A*100, final_pfr$X_A*100))
cat(sprintf("║ Yield Y_P (%%)                                ║  %5.1f  ║   %5.1f  ║ %4.1f ║\n",
    final_batch$Y_P*100, cstr_ss$Y_P*100, final_pfr$Y_P*100))
cat(sprintf("║ Selectivity S_P (%%)                          ║  %5.2f  ║   %5.2f  ║ %4.2f ║\n",
    final_batch$S_P*100, cstr_ss$S_P*100, final_pfr$S_P*100))
cat(sprintf("║ Product Purity (mol%%)                        ║  %5.2f  ║   %5.2f  ║ %4.2f ║\n",
    final_batch$purity*100, cstr_ss$purity*100, final_pfr$purity*100))
cat(sprintf("║ Paracetamol concentration (mol/L)            ║  %5.4f ║   %6.4f ║ %4.4f ║\n",
    final_batch$CP, cstr_ss$CP, final_pfr$CP))
cat(sprintf("║ Impurity concentration (mol/L)               ║  %5.4f ║   %6.4f ║ %4.4f ║\n",
    final_batch$CI, cstr_ss$CI, final_pfr$CI))
cat(sprintf("║ Production rate (mol/min)                  ║  %6.2f  ║   %6.2f  ║  %6.2f  ║\n",
            final_batch$productivity, final_cstr$production_rate, final_pfr$production_rate))
cat("╚══════════════════════════════════════════════╩═════════╩══════════╩══════╝\n")

# ===========================================================================
# WRITE RESULTS TO CSV
# ===========================================================================
results_summary <- data.frame(
  Metric = c("Temperature_K", "Temperature_C",
             "ReactionTime_or_SpaceTime_min",
             "Conversion_pct", "Yield_pct", "Selectivity_pct",
             "Purity_mol_pct", "Paracetamol_mol_L", "Impurity_mol_L"),
  Batch = c(params$T_batch, params$T_batch-273, params$t_batch_end,
            round(final_batch$X_A*100,2), round(final_batch$Y_P*100,2),
            round(final_batch$S_P*100,3), round(final_batch$purity*100,3),
            round(final_batch$CP,4), round(final_batch$CI,5)),
  CSTR  = c(params$T_cstr, params$T_cstr-273, params$tau_cstr,
            round(cstr_ss$X_A*100,2), round(cstr_ss$Y_P*100,2),
            round(cstr_ss$S_P*100,3), round(cstr_ss$purity*100,3),
            round(cstr_ss$CP,4), round(cstr_ss$CI,5)),
  PFR   = c(params$T_pfr, params$T_pfr-273, round(tau_pfr,2),
            round(final_pfr$X_A*100,2), round(final_pfr$Y_P*100,2),
            round(final_pfr$S_P*100,3), round(final_pfr$purity*100,3),
            round(final_pfr$CP,4), round(final_pfr$CI,5))
)
write.csv(results_summary, "Output/reactor_comparison_summary.csv", row.names = FALSE)

# Write full profiles
write.csv(batch_results, "Output/batch_profile.csv", row.names = FALSE)
write.csv(cstr_dynamic,  "Output/cstr_dynamic_profile.csv", row.names = FALSE)
write.csv(pfr_results,   "Output/pfr_profile.csv", row.names = FALSE)
write.csv(cstr_sweep,    "Output/cstr_tau_sweep.csv", row.names = FALSE)
write.csv(pfr_sweep,     "Output/pfr_volume_sweep.csv", row.names = FALSE)

cat("\n✓ CSV data files written to output/\n")

# ===========================================================================
# KEY ENGINEERING CONCLUSIONS
# ===========================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                  KEY ENGINEERING CONCLUSIONS                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  1. PFR achieves comparable conversion and yield for the       ║\n")
cat("║     same volume / space time — ideal plug flow suppresses       ║\n")
cat("║     back-mixing, making it most efficient.                      ║\n")
cat("║                                                                  ║\n")
cat("║  2. CSTR gives the lowest conversion at equivalent τ due to     ║\n")
cat("║     back-mixing (operates at exit concentration throughout).    ║\n")
cat("║     Needs larger τ to match PFR conversion.                     ║\n")
cat("║                                                                  ║\n")
cat("║  3. Batch reactor matches PFR kinetically (identical ODE form)  ║\n")
cat("║     and has higher conversion, suffers from non-productive       ║\n")
cat("║           turn around time and limited throughput                ║\n")
cat("║              suited for small-scale production.                 ║\n")
cat("║                                                                  ║\n")
cat("║  4. Purity is highest in the Batch/PFR since lower T allows     ║\n")
cat("║     better selectivity (side rxn has higher Ea → penalised     ║\n")
cat("║     more at high T). CSTR at 60°C sacrifices some purity for    ║\n")
cat("║     throughput.                                                  ║\n")
cat("║                                                                  ║\n")
cat("║  5. For industrial scale: PFR (or tube reactor) preferred;      ║\n")
cat("║     for flexible multi-product: batch is favoured.              ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

cat("\n✓ Study complete.\n")

