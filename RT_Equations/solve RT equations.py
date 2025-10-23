"""
Simple 0D Rothwarf-Taylor Equations Solver
===========================================
Models quasiparticle and phonon dynamics in supercooled aluminum following
an instantaneous phonon burst (e.g., from cosmic ray muon event).

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
import pickle
import sys

# Add parent directory to path for importing constants
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from constants import (
    eV_to_J, k_B_eV_K,
    Delta_Al_eV, N_0_per_eV, N_cp,
    d_Al, A_qubit, V_Al,
    T_bath, R, B, tau_esc, tau_anh, tau_trap, f_pb,
    n_qp_eq, x_qp_eq, n_ph_eq,
    t_max_RT, n_points_RT,
    file_phonon_injection, dir_rt_equations
)

#==============================================================================
# LOAD PHONON INJECTION DATA
#==============================================================================
print("=" * 80)
print("LOADING PHONON INJECTION DATA")
print("=" * 80)

phonon_data_file = file_phonon_injection
try:
    with open(phonon_data_file, 'rb') as f:
        phonon_data = pickle.load(f)

    phonon_injection_func = phonon_data['phonon_injection_func']
    phonon_metadata = phonon_data['metadata']

    print(f"Loaded phonon injection function from: {phonon_data_file}")
    print(f"  Phonon arrival time: {phonon_data['t_arrival']*1e9:.2f} ns")
    print(f"  Peak phonon flux: {phonon_data['peak_flux']:.2e} m^-3 s^-1")
    print(f"  Total phonon density: {phonon_data['total_phonon_density']:.2e} m^-3")
    print(f"  Energy received at qubit: {phonon_data['energy_received']/1.602e-19:.2e} eV")

    USE_PHONON_INJECTION = True

except FileNotFoundError:
    print(f"X Phonon data file not found: {phonon_data_file}")
    print("  Run 'Phonon Defusion Solver.py' first to generate phonon injection data")
    print("  Continuing with original static burst model...")
    USE_PHONON_INJECTION = False
    phonon_injection_func = None

#==============================================================================
# PHYSICAL CONSTANTS AND PARAMETERS
#==============================================================================

d_film = d_Al
A_film = A_qubit
V_film = V_Al
Delta_Al = Delta_Al_eV
k_B = k_B_eV_K

if USE_PHONON_INJECTION and phonon_injection_func is not None:
    d_film = phonon_metadata['Al_film_thickness']
    A_film = phonon_metadata['qubit_area']
    V_film = phonon_metadata['Al_volume']
    print(f"\nUpdated geometry from phonon data:")
    print(f"  Film thickness: {d_film*1e9:.1f} nm")
    print(f"  Film area: {A_film*1e12:.1f} um^2")
    print(f"  Film volume: {V_film*1e18:.2f} um^3")

#==============================================================================
# PHONON BURST PARAMETERS (for fallback static model)
#==============================================================================
try:
    from constants import file_energy_deposit
    with open(file_energy_deposit, 'rb') as f:
        energy_deposit_data = pickle.load(f)
    E_deposited_Si = energy_deposit_data['energy_deposited_joules']
    muon_energy_from_bethe = energy_deposit_data['muon_energy_mev']
except:
    E_deposited_Si = 1e6 * eV_to_J
    muon_energy_from_bethe = 1.0

eta_transmission = 0.01
E_transmitted = E_deposited_Si * eta_transmission

E_phonon = 2 * Delta_Al * eV_to_J
N_phonons = E_transmitted / E_phonon
n_ph_burst = N_phonons / V_film

#==============================================================================
# INITIAL CONDITIONS
#==============================================================================
if USE_PHONON_INJECTION:
    n_qp_0 = n_qp_eq
    n_ph_0 = n_ph_eq
    print("\nUSING TIME-DEPENDENT PHONON INJECTION")
else:
    N_qp_created = 2 * N_phonons
    n_qp_from_burst = N_qp_created / V_film
    n_qp_0 = n_qp_eq + n_qp_from_burst
    n_ph_0 = n_ph_eq
    print("\nUSING ORIGINAL STATIC BURST MODEL")

print("=" * 80)
print("0D ROTHWARF-TAYLOR EQUATIONS SOLVER")
print("=" * 80)
print("\nSYSTEM PARAMETERS:")
print(f"  Temperature: {T_bath * 1e3:.1f} mK")
print(f"  Al gap: {Delta_Al * 1e6:.1f} ueV")
print(f"  Film volume: {V_film * 1e18:.2f} um^3")
print(f"  Cooper pair density: {N_cp:.2e} m^-3")

print("\nRATE PARAMETERS:")
print(f"  Recombination R: {R:.2e} m^3/s")
print(f"  Pair-breaking B: {B:.2e} m^3/s")
print(f"  Phonon escape time: {tau_esc * 1e9:.2f} ns")
print(f"  Anharmonic decay time: {tau_anh * 1e9:.2f} ns")
print(f"  QP trapping time: {tau_trap * 1e6:.2f} us")
print(f"  PB phonon fraction: {f_pb:.4f}")

print("\nEQUILIBRIUM STATE:")
print(f"  Background n_qp: {n_qp_eq:.2e} m^-3")
print(f"  Background x_qp: {x_qp_eq:.2e}")
print(f"  Background n_ph: {n_ph_eq:.2e} m^-3")

if not USE_PHONON_INJECTION:
    print("\nPHONON BURST & PAIR-BREAKING:")
    print(f"  Muon energy (from Bethe-Bloch): {muon_energy_from_bethe:.1f} MeV ({muon_energy_from_bethe/1000:.2f} GeV)")
    print(f"  Energy deposited in Si: {E_deposited_Si / eV_to_J / 1e6:.2f} MeV")
    print(f"  Transmitted energy to Al: {E_transmitted / eV_to_J / 1e3:.2f} keV")
    print(f"  Phonons arrived: {N_phonons:.2e}")
    print(f"  Cooper pairs broken: {N_phonons:.2e}")
    print(f"  Quasiparticles created: {N_qp_created:.2e}")
    print(f"  QP burst density: {n_qp_from_burst:.2e} m^-3")

print("\nINITIAL CONDITIONS:")
print(f"  n_qp(0): {n_qp_0:.2e} m^-3")
print(f"  x_qp(0): {n_qp_0 / (2 * N_0_per_eV * Delta_Al):.2e}")
print(f"  n_ph(0): {n_ph_0:.2e} m^-3")

#==============================================================================
# ROTHWARF-TAYLOR EQUATIONS
#==============================================================================
def rothwarf_taylor(t, y):
    """
    0D Rothwarf-Taylor equations with trapping and phonon injection

    y[0] = n_qp: quasiparticle density [m⁻³]
    y[1] = n_ph: phonon density (E > 2Δ) [m⁻³]

    Returns: [dn_qp/dt, dn_ph/dt]
    """
    n_qp = max(y[0], 0)
    n_ph = max(y[1], 0)

    # Quasiparticle evolution
    # Factors of 2 account for two quasiparticles per Cooper pair
    dn_qp_dt = (
        -2 * R * n_qp**2
        + 2 * B * n_ph
        - (n_qp - n_qp_eq) / tau_trap
    )

    # Phonon evolution (E > 2Δ only)
    dn_ph_dt = (
        + f_pb * (R / 2) * (2 * n_qp**2)
        - B * n_ph
        - n_ph / tau_esc
        - n_ph / tau_anh
    )

    if USE_PHONON_INJECTION and phonon_injection_func is not None:
        phonon_source = phonon_injection_func(t)
        dn_ph_dt += phonon_source

    return [dn_qp_dt, dn_ph_dt]

#==============================================================================
# SOLVE THE EQUATIONS
#==============================================================================
t_start = 0
t_end = t_max_RT
y0 = [n_qp_0, n_ph_0]

print("\n" + "=" * 80)
print("SOLVING COUPLED EQUATIONS...")
print("=" * 80)

# Solve using stiff ODE solver
solution = solve_ivp(
    rothwarf_taylor,
    t_span=(t_start, t_end),
    y0=y0,
    method='BDF',
    dense_output=True,
    rtol=1e-8,
    atol=1e-10
)

t_eval = np.logspace(-11, np.log10(t_end), n_points_RT)
y_eval = solution.sol(t_eval)
n_qp = y_eval[0]
n_ph = y_eval[1]

if solution.success:
    print("Integration successful")
else:
    print(f"Warning: {solution.message}")

#==============================================================================
# CALCULATE STATISTICS
#==============================================================================
print("\n" + "=" * 80)
print("RESULTS")
print("=" * 80)

# Phonon statistics
idx_ph_max = np.argmax(n_ph)
n_ph_max = n_ph[idx_ph_max]
t_ph_max = t_eval[idx_ph_max]

if n_ph_max > 0:
    try:
        idx_ph_decay = np.where(n_ph[idx_ph_max:] < n_ph_max / np.e)[0]
        t_ph_1e = t_eval[idx_ph_max + idx_ph_decay[0]]
    except IndexError:
        t_ph_1e = np.nan
else:
    t_ph_1e = np.nan

print("\nPHONON POPULATION:")
print(f"  Peak density: {n_ph_max:.2e} m^-3")
print(f"  Peak time: {t_ph_max * 1e9:.2f} ns")
print(f"  1/e decay time: {t_ph_1e * 1e9:.2f} ns")
print(f"  Final density: {n_ph[-1]:.2e} m^-3")

# Quasiparticle statistics
n_qp_peak = n_qp_0
x_qp_peak = n_qp_peak / (2 * N_0_per_eV * Delta_Al)

try:
    idx_qp_relaxed = np.where(n_qp < 2 * n_qp_eq)[0]
    t_qp_relax = t_eval[idx_qp_relaxed[0]]
except IndexError:
    t_qp_relax = np.nan

excess_qp_initial = n_qp_0 - n_qp_eq
target_qp = n_qp_eq + excess_qp_initial / np.e
try:
    idx_qp_1e = np.where(n_qp < target_qp)[0]
    t_qp_1e = t_eval[idx_qp_1e[0]]
except IndexError:
    t_qp_1e = np.nan

print("\nQUASIPARTICLE POPULATION:")
print(f"  Equilibrium density: {n_qp_eq:.2e} m^-3")
print(f"  Equilibrium x_qp: {x_qp_eq:.2e}")
print(f"  Burst density (t=0): {n_qp_peak:.2e} m^-3")
print(f"  Burst x_qp: {x_qp_peak:.2e}")
print(f"  1/e decay time: {t_qp_1e * 1e6:.2f} us")
print(f"  Final density: {n_qp[-1]:.2e} m^-3")
print(f"  Final x_qp: {n_qp[-1] / (2*N_0_per_eV*Delta_Al):.2e}")

# Total particles
N_qp_peak = n_qp_peak * V_film
N_ph_peak = n_ph_max * V_film

print("\nTOTAL PARTICLE NUMBERS:")
print(f"  Initial quasiparticles (burst): {N_qp_peak:.2e}")
print(f"  Peak phonons: {N_ph_peak:.2e}")

print("\n" + "=" * 80)

#==============================================================================
# PLOTTING
#==============================================================================
output_dir = dir_rt_equations
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if USE_PHONON_INJECTION:
    plot_path = os.path.join(output_dir, 'RT_0D_with_phonon_injection.png')
else:
    plot_path = os.path.join(output_dir, 'RT_0D_static_burst.png')


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=False)

# Plot 1: Phonon population
ax1.loglog(t_eval * 1e9, n_ph, 'b-', linewidth=2.5, label='Phonon density n_ph')
ax1.axhline(n_ph_eq, color='gray', linestyle='--', linewidth=1.5,
            label='Background level', alpha=0.7)
if not np.isnan(t_ph_1e):
    ax1.axvline(t_ph_1e * 1e9, color='r', linestyle='--', linewidth=1.5,
                alpha=0.5, label=f'1/e decay time = {t_ph_1e*1e9:.1f} ns')
    ax1.plot(t_ph_max * 1e9, n_ph_max, 'ro', markersize=10,
             label=f'Peak = {n_ph_max:.2e} m⁻³', zorder=5)

ax1.set_xlabel('Time (ns)', fontsize=13, fontweight='bold')
ax1.set_ylabel('Phonon Density n_ph (m⁻³)', fontsize=13, fontweight='bold')
title_suffix = "with Phonon Injection" if USE_PHONON_INJECTION else "Static Burst Model"
ax1.set_title(f'Phonon Population Dynamics (E > 2Δ)\n0D Rothwarf-Taylor Model {title_suffix}',
              fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3, which='both')
ax1.legend(loc='best', fontsize=11)
ax1.set_xlim(left=max(t_eval[0] * 1e9, 1e-3))

# Plot 2: Quasiparticle population
ax2.semilogx(t_eval * 1e6, n_qp, 'r-', linewidth=2.5, label='Quasiparticle density n_qp')
ax2.axhline(n_qp_eq, color='gray', linestyle='--', linewidth=1.5,
            label=f'Equilibrium = {n_qp_eq:.2e} m⁻³', alpha=0.7)
if not np.isnan(t_qp_1e):
    ax2.axvline(t_qp_1e * 1e6, color='b', linestyle='--', linewidth=1.5,
                alpha=0.5, label=f'1/e decay time = {t_qp_1e*1e6:.2f} μs')
ax2.plot(t_eval[0] * 1e6, n_qp_peak, 'go', markersize=10,
         label=f'Initial burst = {n_qp_peak:.2e} m⁻³', zorder=5)

ax2.set_xlabel('Time (μs)', fontsize=13, fontweight='bold')
ax2.set_ylabel('Quasiparticle Density n_qp (m⁻³)', fontsize=13, fontweight='bold')
ax2.set_title(f'Quasiparticle Population Dynamics\n0D Rothwarf-Taylor Model {title_suffix}',
              fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, which='both')
ax2.legend(loc='best', fontsize=11)
ax2.set_xlim(left=max(t_eval[0] * 1e6, 1e-6))

# Plot 3: Normalized quasiparticle density
x_qp = n_qp / (2 * N_0_per_eV * Delta_Al)
ax3.semilogx(t_eval * 1e6, x_qp, 'purple', linewidth=2.5, label='Normalized x_qp')
ax3.axhline(x_qp_eq, color='gray', linestyle='--', linewidth=1.5,
            label=f'Equilibrium x_qp = {x_qp_eq:.2e}', alpha=0.7)
if not np.isnan(t_qp_1e):
    ax3.axvline(t_qp_1e * 1e6, color='b', linestyle='--', linewidth=1.5,
                alpha=0.5, label=f'1/e decay time = {t_qp_1e*1e6:.2f} μs')
ax3.plot(t_eval[0] * 1e6, x_qp[0], 'go', markersize=10,
         label=f'Initial x_qp = {x_qp[0]:.2e}', zorder=5)
ax3.plot(t_eval[-1] * 1e6, x_qp[-1], 'ro', markersize=10,
         label=f'Final x_qp = {x_qp[-1]:.2e}', zorder=5)

ax3.set_xlabel('Time (μs)', fontsize=13, fontweight='bold')
ax3.set_ylabel('Normalized QP Density x_qp', fontsize=13, fontweight='bold')
ax3.set_title(f'Normalized Quasiparticle Density (x_qp = n_qp / 2N₀Δ)\n0D Rothwarf-Taylor Model {title_suffix}',
              fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3, which='both')
ax3.legend(loc='best', fontsize=11)
ax3.set_xlim(left=max(t_eval[0] * 1e6, 1e-6))

plt.tight_layout(pad=3.0)
plt.savefig(plot_path, dpi=300, bbox_inches='tight')

print(f"\nPlot saved to: {plot_path}")

# ============================================================================
# EXPORT TIME-DEPENDENT x_qp(t) FOR T1 CALCULATIONS
# ============================================================================
print("\n" + "=" * 80)
print("EXPORTING x_qp(t) DATA FOR T1 CALCULATIONS")
print("=" * 80)

x_qp_t = n_qp / (2 * N_0_per_eV * Delta_Al)

xqp_data = {
    't_array': t_eval,
    'x_qp': x_qp_t,
    'n_qp': n_qp,
    'n_ph': n_ph,
    'x_qp_eq': x_qp_eq,
    'n_qp_eq': n_qp_eq,
    'metadata': {
        'Delta_Al': Delta_Al,
        'T_bath': T_bath,
        'N_0_per_eV': N_0_per_eV,
        'V_film': V_film,
        'd_film': d_film,
        'A_film': A_film,
        'use_phonon_injection': USE_PHONON_INJECTION,
    }
}

output_file = f'{dir_rt_equations}/xqp_vs_time_data.pkl'
with open(output_file, 'wb') as f:
    pickle.dump(xqp_data, f)

print(f"\nTime-dependent x_qp(t) data saved to: {output_file}")
print(f"\nData properties:")
print(f"  Time range: {t_eval[0]*1e9:.2f} ns to {t_eval[-1]*1e3:.2f} ms")
print(f"  Number of time points: {len(t_eval)}")
print(f"  Initial x_qp: {x_qp_t[0]:.2e}")
print(f"  Final x_qp: {x_qp_t[-1]:.2e}")
print(f"  Peak x_qp: {np.max(x_qp_t):.2e}")
print(f"  Equilibrium x_qp: {x_qp_eq:.2e}")

print("\n" + "=" * 80)
print("SIMULATION COMPLETE")
print("=" * 80)
