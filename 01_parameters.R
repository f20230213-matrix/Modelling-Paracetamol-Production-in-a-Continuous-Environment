# =============================================================================
# Paracetamol (Acetaminophen) Synthesis - Reaction Parameters
# =============================================================================
# Reaction: 4-aminophenol + Acetic anhydride → Paracetamol + Acetic acid
#
# This is a liquid-phase acetylation reaction modelled as:
#   A + B → P + S
# where:
#   A = 4-aminophenol  (limiting reagent)
#   B = Acetic anhydride
#   P = Paracetamol
#   S = Acetic acid (byproduct)
#
# A parallel side reaction produces an N,O-diacetyl impurity:
#   A + 2B → I + 2S   (undesired)
# =============================================================================

params <- list(

  # --- Kinetic constants (Arrhenius) ---
  # Main reaction:  A + B -> P + S
  k1_ref  = 0.08,          # L/(mol·min) at T_ref
  Ea1     = 55000,          # J/mol  (activation energy, main rxn)

  # Side reaction: A + B -> I (diacetyl impurity)
  k2_ref  = 0.005,          # L/(mol·min) at T_ref
  Ea2     = 70000,          # J/mol  (activation energy, side rxn - higher Ea => more sensitive to T)

  T_ref   = 300,            # K  (reference temperature = 27°C)
  R_gas   = 8.314,          # J/(mol·K)

  # --- Operating temperatures ---
  T_batch = 323,            # K  (50°C)  - moderate to control side rxn
  T_cstr  = 333,            # K  (60°C)  - slightly higher for throughput
  T_pfr   = 328,            # K  (55°C)  - intermediate

  # --- Feed / initial conditions ---
  CA0 = 1.0,    # mol/L  initial 4-aminophenol concentration
  CB0 = 1.2,    # mol/L  initial acetic anhydride (slight excess)
  CP0 = 0.0,
  CI0 = 0.0,    # impurity

  # --- Batch reactor ---
  t_batch_end = 60,         # min  (total reaction time)
  V_batch_total = 100,      # L    (reactor volume)

  # --- CSTR ---
  tau_cstr    = 30,         # min  (residence time)
  V_cstr      = 100,        # L    (reactor volume)

  # --- PFR ---
  V_pfr_total = 100,        # L    (reactor volume)
  Q_pfr       = 3.0,        # L/min (volumetric flow rate)

  # --- Selectivity / Purity target ---
  target_conv = 0.90        # 90% conversion target
)

# Arrhenius rate constant at given T
k_arrhenius <- function(k_ref, Ea, T, T_ref = 300, R = 8.314) {
  k_ref * exp(-Ea / R * (1/T - 1/T_ref))
}

cat("Parameters loaded successfully.\n")
cat(sprintf("  k1(T_batch=%.0fK) = %.5f L/(mol·min)\n",
            params$T_batch, k_arrhenius(params$k1_ref, params$Ea1, params$T_batch)))
cat(sprintf("  k2(T_batch=%.0fK) = %.6f L/(mol·min)\n",
            params$T_batch, k_arrhenius(params$k2_ref, params$Ea2, params$T_batch)))

