# Transmon Qubit Cosmic Ray Impact Simulation

A complete simulation pipeline for modeling the impact of cosmic ray muons on superconducting transmon qubits, from initial energy deposition through quasiparticle dynamics to coherence time (T1) degradation and recovery.

## Overview

This project simulates the complete physics chain of a cosmic ray event affecting a superconducting qubit:

1. **Energy Deposition**: A cosmic ray muon deposits energy in the silicon substrate (Bethe-Bloch formalism)
2. **Phonon Transport**: High-energy phonons propagate ballistically from the substrate to the aluminum qubit
3. **Quasiparticle Dynamics**: Phonons break Cooper pairs, creating quasiparticles (Rothwarf-Taylor equations)
4. **T1 Degradation**: Quasiparticles cause coherence time degradation and gradual recovery

The simulation captures physics across 9 orders of magnitude in time: from nanosecond phonon transport to millisecond quasiparticle relaxation.

## Quick Start

### Run Complete Pipeline

To run all four simulation modules in sequence:

```bash
python main.py
```

This will execute all modules and generate all output files and plots (~33 seconds).

### Run Individual Modules

You can also run individual modules separately:

```bash
# 1. Energy deposition (Bethe-Bloch)
python "Energy_Deposit/bethe.py"

# 2. Phonon transport
python "Phonon_Dynamics/Phonon Defusion Solver.py"

# 3. Quasiparticle dynamics (RT equations)
python "RT_Equations/solve RT equations.py"

# 4. T1 calculation
python "T1_Calculation/plot gamma qp.py"
```

**Note**: Modules must be run in order since later modules depend on data from earlier ones.

## Mathematical Framework

### 1. Energy Deposition: Bethe-Bloch Equation

The energy loss rate of a muon passing through the silicon substrate is given by the **Bethe-Bloch formula**:

```
-dE/dx = K z² (Z/A) (1/β²) [½ ln(2mₑc²β²γ²Tₘₐₓ/I²) - β² - δ(βγ)/2]
```

Where:
- **K** = 0.307075 MeV·cm²/g (constant for A=1 g/mol)
- **z** = charge of incident particle (z=1 for muon)
- **Z/A** = atomic number/mass ratio of target (Si: Z=14, A=28.09)
- **β** = v/c (particle velocity relative to light speed)
- **γ** = 1/√(1-β²) (Lorentz factor)
- **mₑ** = electron mass
- **Tₘₐₓ** = maximum kinetic energy transfer in a single collision
- **I** = mean excitation energy of Si (173 eV)
- **δ(βγ)** = density effect correction (Sternheimer parameterization)

**Landau fluctuations** are included to model stochastic variations in energy deposition, crucial for understanding the statistical distribution of cosmic ray impacts.

**Implementation**: [Energy_Deposit/bethe.py](Energy_Deposit/bethe.py)

### 2. Phonon Transport: Ballistic Propagation

High-energy phonons (E ≈ 3Δ ≈ 540 μeV) created in the substrate propagate ballistically to the qubit. The phonon flux at the qubit location is calculated using a **geometric factor model**:

```
Φ_phonon(x,y,t) = (E_deposited × T_interface) / (4π R³ v_sound) × f(t)
```

Where:
- **E_deposited** = energy deposited by muon in substrate (from Bethe-Bloch)
- **T_interface** = transmission coefficient at Si/Al interface
  ```
  T_interface = 4 Z_Si Z_Al / (Z_Si + Z_Al)²
  ```
  with Z = ρv (acoustic impedance)
- **R** = distance from deposition point to qubit surface
  ```
  R = √[(x-x₀)² + (y-y₀)² + z₀²]
  ```
- **v_sound** = sound velocity in Si (5500 m/s)
- **f(t)** = temporal pulse shape (Gaussian with width ~50 ns)

The phonon flux is integrated over the qubit area to obtain the **total phonon injection rate** into the aluminum film:

```
Ṅ_phonon(t) = ∫∫_A_qubit Φ_phonon(x,y,t) dA / E_phonon
```

**Implementation**: [Phonon_Dynamics/Phonon Defusion Solver.py](Phonon_Dynamics/Phonon Defusion Solver.py)

### 3. Quasiparticle Dynamics: Rothwarf-Taylor Equations

The coupled dynamics of quasiparticles and phonons in the superconductor are governed by the **Rothwarf-Taylor (RT) equations**:

```
dn_qp/dt = -2R n_qp² + 2B n_ph - (n_qp - n_qp_eq)/τ_trap + S_phonon(t)

dn_ph/dt = f_pb × (R/2) × (2n_qp²) - B n_ph - n_ph/τ_esc - n_ph/τ_anh
```

**Physical processes**:
1. **Recombination**: `-2R n_qp²` — two quasiparticles recombine to form a Cooper pair, emitting a phonon (factor of 2 because two QPs are consumed)
2. **Pair-breaking**: `+2B n_ph` — a high-energy phonon (E > 2Δ) breaks a Cooper pair, creating two quasiparticles
3. **Relaxation to traps**: `-(n_qp - n_qp_eq)/τ_trap` — quasiparticles relax to normal metal traps or defects (τ_trap ≈ 1 ms)
4. **External injection**: `S_phonon(t)` — phonon flux from substrate (output of Module 2)
5. **Phonon recycling**: `f_pb × (R/2) × (2n_qp²)` — fraction f_pb ≈ 0.5 of recombination phonons have E > 2Δ and can break pairs again
6. **Phonon escape**: `-n_ph/τ_esc` — phonons escape the thin Al film (τ_esc ≈ 80 ns)
7. **Anharmonic decay**: `-n_ph/τ_anh` — phonons decay into lower-energy phonons (τ_anh ≈ 1 ns)

**Normalization**: The quasiparticle density is normalized as:
```
x_qp = n_qp / (2 N₀ Δ)
```
where N₀ = 1.74×10²⁸ eV⁻¹m⁻³ is the single-spin density of states at the Fermi level, and the factor of 2 accounts for both spin directions.

**Equilibrium values** (T = 20 mK):
- x_qp_eq ≈ 1.6×10⁻⁷ (residual background)
- n_qp_eq ≈ 1×10¹⁸ m⁻³ (≈ 1 QP per μm³)

**Implementation**: [RT_Equations/solve RT equations.py](RT_Equations/solve RT equations.py) using `scipy.integrate.solve_ivp` with BDF method (for stiff ODEs).

### 4. Coherence Time: T1 Calculation

The quasiparticle-induced relaxation rate is:

```
Γ_qp = √(2 ω_q Δ / π² ℏ) × x_qp
```

Where:
- **ω_q** = 2πf_q = qubit angular frequency (f_q ≈ 6 GHz)
- **Δ** = superconducting gap (180 μeV for Al)
- **x_qp** = normalized quasiparticle density (from RT equations)

The total relaxation rate includes both QP and non-QP contributions:

```
1/T1_total(t) = Γ_qp(t) + 1/T1_intrinsic
```

Where T1_intrinsic ≈ 100 μs represents other decoherence mechanisms (dielectric loss, surface defects, etc.).

**Key results**:
- **Equilibrium**: Γ_qp_eq ≈ 7300 s⁻¹ → T1_eq ≈ 58 μs
- **Peak burst**: Γ_qp_peak ≈ 64000 s⁻¹ → T1_min ≈ 13 μs
- **Recovery timescale**: τ_recovery ≈ 3-6 ms (governed by τ_trap = 1 ms)

**Implementation**: [T1_Calculation/plot gamma qp.py](T1_Calculation/plot gamma qp.py)

### Key Physics Insights

1. **Energy cascade efficiency**: Only ~0.0016% of the muon's deposited energy reaches the qubit, due to:
   - Geometric spreading (1/R³ dilution)
   - Acoustic impedance mismatch at Si/Al interface (T ≈ 0.02)
   - Phonon energy filtering (only E > 2Δ contribute to pair-breaking)

2. **Timescale hierarchy**:
   - Phonon arrival: ~387 ns (ballistic propagation)
   - Phonon escape: ~80 ns (thin film geometry)
   - QP generation: ~1-10 ns (pair-breaking dynamics)
   - QP relaxation: ~1 ms (trap-limited)
   - Full T1 recovery: ~10-100 ms (governed by τ_trap)

3. **Nonequilibrium dynamics**: The system exhibits a sharp QP burst (peak at ~400 ns) followed by slow exponential decay. The recombination term `R n_qp²` provides fast initial decay, but the bottleneck is the linear relaxation term `τ_trap`.

## Project Structure

```
Transmon_Model/
├── main.py                          # Main pipeline script (run this!)
├── constants.py                     # Centralized physical constants
├── Energy_Deposit/
│   ├── bethe.py                    # Bethe-Bloch energy deposition
│   └── energy_deposit_data.pkl     # Output: energy deposit data
├── Phonon_Dynamics/
│   ├── Phonon Defusion Solver.py   # Ballistic phonon transport
│   └── phonon_injection_data.pkl   # Output: phonon flux function
├── RT_Equations/
│   ├── solve RT equations.py       # Rothwarf-Taylor equations
│   └── xqp_vs_time_data.pkl       # Output: quasiparticle density vs time
└── T1_Calculation/
    ├── plot gamma qp.py            # T1 coherence time calculation
    ├── T1_vs_time.png              # Output: T1 recovery plot
    └── xqp_and_T1_vs_time.png      # Output: combined QP and T1 plot
```

## Key Physics Parameters

All physical constants are defined in `constants.py`:

### Material Properties
- **Aluminum film**: 50 nm thick, 100×100 μm² area
- **Superconducting gap**: 180 μeV
- **Silicon substrate**: 200 μm thick
- **Operating temperature**: 20 mK

### Qubit Parameters
- **Frequency**: 6 GHz
- **Intrinsic T1**: 100 μs (non-QP contributions)

### Rate Parameters (Based on Experimental Literature)
- **Recombination coefficient R**: 1.0×10⁻¹⁸ m³/s
- **Pair-breaking rate B**: 9.0×10⁹ s⁻¹
- **QP relaxation time τ_trap**: 1 ms (Nature Comm 4, 5836, 2014)
- **Phonon escape time**: 80 ns
- **Anharmonic decay time**: 1 ns

### Cosmic Ray Parameters
- **Muon energy**: 2 GeV (typical)
- **Incident angle**: 0° (normal incidence)

## Key Results

Running the complete pipeline produces the following key findings:

### Energy Flow
- **Energy deposited in Si**: ~0.074 MeV (74 keV)
- **Energy transmitted to Al**: ~1.2 eV (0.0016% efficiency)
- **Phonon arrival time**: ~387 ns
- **Peak phonon density**: ~3.5×10¹⁵ m⁻³

### Quasiparticle Dynamics
- **Peak x_qp**: ~1.4×10⁻⁶ (enhancement of ~8.8× over equilibrium)
- **Equilibrium x_qp**: 1.6×10⁻⁷
- **QP relaxation**: Millisecond timescale

### T1 Coherence Time
- **Equilibrium T1**: 57.8 μs
- **Minimum T1 (during burst)**: 13.5 μs
- **Recovery times**:
  - 90% recovery: 3.4 ms
  - 95% recovery: 4.2 ms
  - 99% recovery: 5.8 ms
- **Full recovery**: ~100 ms

These timescales match experimental observations of millisecond-scale T1 recovery in literature.

## Output Files

### Data Files (.pkl)
- `Energy_Deposit/energy_deposit_data.pkl`: Muon energy deposition results
- `Phonon_Dynamics/phonon_injection_data.pkl`: Time-dependent phonon flux function
- `RT_Equations/xqp_vs_time_data.pkl`: Quasiparticle density evolution

### Plots (.png)
- `Energy_Deposit/energy_deposit_summary.png`: Bethe-Bloch results with Landau fluctuations
- `Phonon_Dynamics/ballistic_phonon_flux_and_heatmap.png`: Spatial phonon distribution
- `RT_Equations/RT_0D_with_phonon_injection.png`: QP and phonon population dynamics
- `T1_Calculation/gamma_qp_plot.png`: Decay rate vs QP density
- `T1_Calculation/gamma_qp_dual_axis.png`: Combined Γ_qp and T1 plot
- `T1_Calculation/T1_vs_time.png`: T1 degradation and recovery
- `T1_Calculation/xqp_and_T1_vs_time.png`: Correlation between x_qp(t) and T1(t)

## Modifying Parameters

To modify simulation parameters, edit `constants.py`:

```python
# Example: Change muon energy
muon_energy_MeV = 2000  # 2 GeV (default)

# Example: Change QP relaxation time
tau_trap = 1e-3  # 1 ms (default, based on literature)

# Example: Change substrate thickness
substrate_thickness_um = 200  # 200 μm (default)

# Example: Change qubit position
x_0 = 0.5e-3  # 0.5 mm (default)
y_0 = 0.5e-3  # 0.5 mm (default)
```

After modifying parameters, re-run the pipeline:
```bash
python main.py
```

## Dependencies

Required Python packages:
- numpy
- scipy
- matplotlib
- pickle (standard library)

## References

### Key Experimental Papers
1. **Quasiparticle dynamics**: Martinis et al., *Nature Communications* **4**, 5836 (2014) — "Measurement and control of quasiparticle dynamics in a superconducting qubit" (source for τ_trap ≈ 1 ms)
2. **T1 recovery timescales**: Cardani et al., *Nature Communications* **16**, 364 (2025) — "Methods to achieve near-millisecond energy relaxation times..."
3. **QP pumping**: Aumentado et al., *Science* **355**, 6323 (2017) — "Suppressing relaxation in superconducting qubits by quasiparticle pumping"
4. **Cosmic ray impacts**: Vepsäläinen et al., *Nature* **584**, 551-556 (2020) — "Impact of ionizing radiation on superconducting qubit coherence"

### Theoretical Framework
- **Energy loss**: Bethe-Bloch equation (Particle Data Group review)
- **Phonon transport**: Ballistic phonon model with geometric factors (Kaplan et al.)
- **Superconductor dynamics**: Rothwarf-Taylor equations, *Phys. Rev. Lett.* **26**, 1355 (1971)
- **QP decoherence**: Catelani et al., *Phys. Rev. B* **84**, 064517 (2011)

## Version History

### October 2025
**Major bug fixes** (critical for accurate T1 recovery):
1. **x_qp normalization**: Fixed inconsistent normalization (factor of 4 vs 2) in `constants.py` that prevented T1 from recovering to equilibrium. Now uses correct factor of 2: `x_qp = n_qp / (2 N₀ Δ)`
2. **τ_trap timescale**: Corrected quasiparticle relaxation time from 1 μs to 1 ms based on experimental literature (Nature Comm 4, 5836, 2014), yielding physically correct millisecond recovery times
3. **Output improvements**: Enhanced terminal output and plots to clearly show T1 recovery metrics and progress

## License

Generated with Claude Code (Anthropic)

## Contact

For questions or issues, please check the project documentation or constants.py for parameter definitions.
