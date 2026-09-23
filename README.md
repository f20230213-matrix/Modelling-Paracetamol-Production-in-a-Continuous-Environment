# Modelling Paracetamol Production in a Continuous Environment

A reaction-engineering simulation study comparing **batch, continuous stirred-tank reactor (CSTR), and plug-flow reactor (PFR)** configurations for the synthesis of paracetamol.

The project investigates how reactor configuration and operating conditions affect **conversion, yield, selectivity, purity, and productivity**.

> **Note:** This is an academic/process-modelling study based on idealized reactor models. The results are intended for reactor comparison and engineering analysis, not as a plant-ready process design.

---

## Project Overview

Paracetamol (acetaminophen) production has traditionally been associated with batch processing. Continuous manufacturing offers potential advantages such as reduced residence time, steady-state operation, improved throughput, and greater process consistency.

This project compares three reactor configurations:

- **Batch Reactor**
- **Continuous Stirred-Tank Reactor (CSTR)**
- **Plug Flow Reactor (PFR)**

The objective is to understand the effect of reactor configuration on paracetamol synthesis and evaluate the trade-offs between conversion, selectivity, purity, and productivity.

---

## Reaction System

The primary reaction is the acetylation of 4-aminophenol:

$$
\text{4-Aminophenol} + \text{Acetic Anhydride}
\rightarrow
\text{Paracetamol} + \text{Acetic Acid}
$$

A parallel undesired reaction leading to an N,O-diacetyl impurity is also incorporated into the model.

The reaction rates are represented as:

$$
r_1 = k_1 C_A C_B
$$

$$
r_2 = k_2 C_A C_B
$$

where:

- $A$ = 4-aminophenol
- $B$ = acetic anhydride
- $P$ = paracetamol
- $I$ = undesired impurity

The temperature dependence of the rate constants is represented using the Arrhenius equation:

$$
k(T)=k_{\mathrm{ref}}
\exp\left[
-\frac{E_a}{R}
\left(
\frac{1}{T}-\frac{1}{T_{\mathrm{ref}}}
\right)
\right]
$$

---

## Reactor Models

### Batch Reactor

The batch reactor is modeled using transient material balances:

$$
\frac{dC_A}{dt}=-(r_1+r_2)
$$

$$
\frac{dC_B}{dt}=-(r_1+r_2)
$$

$$
\frac{dC_P}{dt}=r_1
$$

$$
\frac{dC_I}{dt}=r_2
$$

The differential equations are solved numerically using the `deSolve` package and the LSODA solver.

### CSTR

The CSTR model includes both steady-state and dynamic startup simulations.

Residence time is defined as:

$$
\tau = \frac{V}{Q}
$$

A residence-time sweep is also performed to determine the operating conditions required to approach the target conversion.

### PFR

The PFR model is formulated using space time:

$$
\tau = \frac{V}{Q}
$$

The concentration profiles are obtained by integrating the reaction balances along reactor space time.

A reactor-volume sweep is performed to examine the volume required to achieve the desired conversion.

---

## Model Parameters

The baseline model uses the following conditions:

| Parameter | Batch | CSTR | PFR |
|---|---:|---:|---:|
| Temperature | 50 °C | 60 °C | 55 °C |
| Reactor / Space Time | 60 min | 30 min | 33.33 min |
| Reactor Volume | 100 L | 100 L | 100 L |
| 4-Aminophenol | 1.0 mol/L | 1.0 mol/L | 1.0 mol/L |
| Acetic Anhydride | 1.2 mol/L | 1.2 mol/L | 1.2 mol/L |

### Kinetic Parameters

| Parameter | Value |
|---|---:|
| $k_{1,\mathrm{ref}}$ | 0.08 L/(mol·min) |
| $E_{a,1}$ | 55 kJ/mol |
| $k_{2,\mathrm{ref}}$ | 0.005 L/(mol·min) |
| $E_{a,2}$ | 70 kJ/mol |
| Reference Temperature | 300 K |
| Target Conversion | 90% |

---

## Results

The baseline simulation gives the following reactor comparison:

| Metric | Batch | CSTR | PFR |
|---|---:|---:|---:|
| Temperature | 50 °C | 60 °C | 55 °C |
| Reaction / Space Time | 60 min | 30 min | 33.33 min |
| Conversion | 99.89% | 88.29% | 99.65% |
| Paracetamol Yield | 91.15% | 79.29% | 90.22% |
| Selectivity | 91.25% | 89.81% | 90.54% |
| Product Purity | 91.25% | 89.81% | 90.54% |
| Paracetamol Concentration | 0.9115 mol/L | 0.7929 mol/L | 0.9022 mol/L |
| Impurity Concentration | 0.0874 mol/L | 0.0899 mol/L | 0.0942 mol/L |

### Conversion Target

The model uses a **90% conversion target**.

The sensitivity analyses indicate approximately:

- **Batch reactor:** ~11 min to reach 90% conversion
- **CSTR:** ~40 min residence time to approach 90% conversion
- **PFR:** ~8.1 min space time to reach 90% conversion

For the modeled flow rate, the PFR conversion target corresponds to approximately **24.3 L reactor volume**.

---

## Analysis Included

The project performs several analyses.

### 1. Concentration Profiles

Concentration profiles are generated for:

- 4-Aminophenol
- Acetic Anhydride
- Paracetamol
- Undesired Impurity

### 2. Reactor Performance Comparison

The reactor configurations are compared using:

- Conversion
- Yield
- Selectivity
- Product purity
- Production rate

### 3. Residence-Time and Volume Sensitivity

The model investigates:

- CSTR conversion vs. residence time
- PFR conversion vs. reactor volume

### 4. Temperature Sensitivity

The effect of temperature on the desired and undesired reaction pathways is investigated using the Arrhenius relationship.

### 5. Summary Comparison

The generated plots provide a visual comparison of reactor performance across the key metrics.

---

## Repository Structure

```text
Modelling-Paracetamol-Production-in-a-Continuous-Environment/
│
├── 01_parameters.R
├── main.R
│
├── R code/
│   ├── 02_batch_reactor.R
│   ├── 03_cstr_model.R
│   ├── 04_pfr_model.R
│   └── 05_comparison_plots.R
│
├── Output/
│   ├── reactor_comparison_summary.csv
│   ├── batch_profile.csv
│   ├── cstr_dynamic_profile.csv
│   ├── cstr_tau_sweep.csv
│   ├── pfr_profile.csv
│   └── pfr_volume_sweep.csv
│
└── Plots/
    ├── 01_concentration_profiles.png
    ├── 02_performance_metrics.png
    ├── 03_sensitivity_analysis.png
    ├── 04_temperature_effects.png
    └── 05_summary_comparison.png
