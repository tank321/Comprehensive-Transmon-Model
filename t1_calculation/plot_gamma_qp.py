"""
Plot of Quasiparticle Decay Rate (Gamma_qp) vs Normalized QP Density (x_qp)
Based on Equation 31 from the model document
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import sys

# Add parent directory to path for importing constants
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from constants import (
    Delta_Al_eV,
    omega_q, f_q,
    T1_intrinsic, Gamma_intrinsic,
    Gamma_qp_prefactor,
    x_qp_typical_residual, x_qp_burst_low, x_qp_burst_high,
    dir_t1_calculation
)

#==============================================================================
# LOAD TIME-DEPENDENT x_qp(t) DATA FROM RT EQUATIONS
#==============================================================================
print("=" * 80)
print("LOADING TIME-DEPENDENT x_qp(t) DATA")
print("=" * 80)

xqp_data_file = 'rt_equations/xqp_vs_time_data.pkl'
try:
    with open(xqp_data_file, 'rb') as f:
        xqp_data = pickle.load(f)

    t_array = xqp_data['t_array']           # Time array [s]
    x_qp_t = xqp_data['x_qp']               # Normalized QP density vs time
    n_qp_t = xqp_data['n_qp']               # Absolute QP density [m^-3]
    x_qp_eq = xqp_data['x_qp_eq']           # Equilibrium x_qp
    metadata = xqp_data['metadata']

    print(f"Loaded x_qp(t) data from: {xqp_data_file}")
    print(f"  Time range: {t_array[0]*1e9:.2f} ns to {t_array[-1]*1e3:.2f} ms")
    print(f"  Number of points: {len(t_array)}")
    print(f"  x_qp range: {np.min(x_qp_t):.2e} to {np.max(x_qp_t):.2e}")
    print(f"  Using phonon injection: {metadata['use_phonon_injection']}")

    USE_TIME_DEPENDENT_DATA = True

except FileNotFoundError:
    print(f"X Time-dependent data file not found: {xqp_data_file}")
    print("  Run 'solve RT equations.py' first to generate x_qp(t) data")
    print("  Continuing with static analysis only...")
    USE_TIME_DEPENDENT_DATA = False

# Calculate prefactor (use imported constants)
prefactor = Gamma_qp_prefactor

print("=" * 70)
print("QUASIPARTICLE DECAY RATE CALCULATION")
print("=" * 70)
print(f"\nTransmon Parameters:")
print(f"  Qubit frequency: {f_q / 1e9:.2f} GHz")
print(f"  Angular frequency omega_q: {omega_q / (2*np.pi*1e9):.2f} x 2pi GHz")
print(f"  Superconducting gap Delta: {Delta_Al_eV * 1e6:.1f} ueV")
print(f"\nPrefactor sqrt(2*omega_q*Delta/(pi^2*hbar)): {prefactor:.2e} s^-1")

x_qp = np.logspace(-10, -3, 1000)

Gamma_qp = prefactor * x_qp
Gamma_total = Gamma_qp + Gamma_intrinsic

Gamma_qp_kHz = Gamma_qp / 1e3
Gamma_total_kHz = Gamma_total / 1e3
T1_qp = 1 / Gamma_qp
T1_total = 1 / Gamma_total
T1_qp_us = T1_qp * 1e6
T1_total_us = T1_total * 1e6

# Use imported typical x_qp values
x_qp_typical = x_qp_typical_residual
# x_qp_burst_low and x_qp_burst_high already imported

Gamma_qp_typical = prefactor * x_qp_typical
Gamma_qp_burst_low = prefactor * x_qp_burst_low
Gamma_qp_burst_high = prefactor * x_qp_burst_high

T1_qp_typical = 1 / Gamma_qp_typical
T1_qp_burst_low = 1 / Gamma_qp_burst_low
T1_qp_burst_high = 1 / Gamma_qp_burst_high

Gamma_total_typical = Gamma_qp_typical + Gamma_intrinsic
Gamma_total_burst_low = Gamma_qp_burst_low + Gamma_intrinsic
Gamma_total_burst_high = Gamma_qp_burst_high + Gamma_intrinsic

T1_total_typical = 1 / Gamma_total_typical
T1_total_burst_low = 1 / Gamma_total_burst_low
T1_total_burst_high = 1 / Gamma_total_burst_high

print(f"\nIntrinsic Loss Mechanisms (non-QP):")
print(f"  T1_intrinsic = {T1_intrinsic * 1e6:.1f} us")
print(f"  Gamma_intrinsic = {Gamma_intrinsic / 1e3:.2f} kHz")

print(f"\nTypical Operating Conditions:")
print(f"  x_qp = {x_qp_typical:.0e} (residual background)")
print(f"    -> Gamma_qp = {Gamma_qp_typical / 1e3:.2f} kHz")
print(f"    -> T1_qp (QP-only) = {T1_qp_typical * 1e6:.1f} us")
print(f"    -> T1_total (QP + intrinsic) = {T1_total_typical * 1e6:.1f} us")

print(f"\n  x_qp = {x_qp_burst_low:.0e} (low-level burst)")
print(f"    -> Gamma_qp = {Gamma_qp_burst_low / 1e3:.2f} kHz")
print(f"    -> T1_qp (QP-only) = {T1_qp_burst_low * 1e6:.1f} us")
print(f"    -> T1_total (QP + intrinsic) = {T1_total_burst_low * 1e6:.1f} us")

print(f"\n  x_qp = {x_qp_burst_high:.0e} (high-level burst)")
print(f"    -> Gamma_qp = {Gamma_qp_burst_high / 1e3:.2f} kHz")
print(f"    -> T1_qp (QP-only) = {T1_qp_burst_high * 1e6:.2f} us")
print(f"    -> T1_total (QP + intrinsic) = {T1_total_burst_high * 1e6:.2f} us")

print("=" * 70)

fig2, ax = plt.subplots(1, 1, figsize=(12, 8))

ax_T1 = ax.twinx()

line1, = ax.loglog(x_qp, Gamma_qp_kHz, 'b-', linewidth=3, label='Gamma_qp (left axis)')

line2a, = ax_T1.loglog(x_qp, T1_qp_us, 'r--', linewidth=2, alpha=0.5, label='T1_qp (QP-only, right axis)')
line2b, = ax_T1.loglog(x_qp, T1_total_us, 'r-', linewidth=3, label='T1_total (right axis)')

# Add intrinsic T1 limit
ax_T1.axhline(T1_intrinsic * 1e6, color='purple', linestyle='--', linewidth=2, alpha=0.7,
              label=f'Intrinsic T1 = {T1_intrinsic * 1e6:.0f} us')

ax.plot(x_qp_typical, Gamma_qp_typical / 1e3, 'go', markersize=14, zorder=5)
ax_T1.plot(x_qp_typical, T1_total_typical * 1e6, 'go', markersize=14, zorder=5)

ax.plot(x_qp_burst_low, Gamma_qp_burst_low / 1e3, 'yo', markersize=14, zorder=5)
ax_T1.plot(x_qp_burst_low, T1_total_burst_low * 1e6, 'yo', markersize=14, zorder=5)

ax.plot(x_qp_burst_high, Gamma_qp_burst_high / 1e3, 'ro', markersize=14, zorder=5)
ax_T1.plot(x_qp_burst_high, T1_total_burst_high * 1e6, 'ro', markersize=14, zorder=5)

ax.annotate(f'Residual\nx_qp = 10-7\nT1 = {T1_total_typical*1e6:.0f} us',
           xy=(x_qp_typical, Gamma_qp_typical / 1e3),
           xytext=(1e-8, 1e1), fontsize=11,
           arrowprops=dict(arrowstyle='->', lw=2, color='green'),
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

ax.annotate(f'Burst\nx_qp = 10-5\nT1 = {T1_total_burst_high*1e6:.1f} us',
           xy=(x_qp_burst_high, Gamma_qp_burst_high / 1e3),
           xytext=(1e-6, 5e2), fontsize=11,
           arrowprops=dict(arrowstyle='->', lw=2, color='red'),
           bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))

ax.axvspan(1e-8, 1e-6, alpha=0.15, color='green', label='Typical residual')
ax.axvspan(1e-6, 1e-4, alpha=0.15, color='orange', label='Burst events')

ax.set_xlabel('Normalized Quasiparticle Density x_qp = n_qp / (2N₀Δ)', fontsize=14, fontweight='bold')
ax.set_ylabel('QP Decay Rate Gamma_qp (kHz)', fontsize=14, fontweight='bold', color='b')
ax_T1.set_ylabel('Coherence Time T1 (us)', fontsize=14, fontweight='bold', color='r')

ax.set_title('Quasiparticle-Induced Decoherence (Eq. 31)\n' +
            r'$\Gamma_{\rm qp} = \sqrt{\frac{2\omega_q\Delta}{\pi^2\hbar}} \cdot x_{\rm qp}$ for 6 GHz transmon, Δ = 180 ueV',
            fontsize=15, fontweight='bold', pad=20)

ax.tick_params(axis='y', labelcolor='b', labelsize=12)
ax_T1.tick_params(axis='y', labelcolor='r', labelsize=12)
ax.tick_params(axis='x', labelsize=12)

ax.grid(True, alpha=0.3, which='both')

# Create combined legend with all elements
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
line_elements = [line1, line2a, line2b]
intrinsic_line = Line2D([0], [0], color='purple', linestyle='--', linewidth=2,
                        label=f'Intrinsic T1 = {T1_intrinsic * 1e6:.0f} us')
patch_elements = [Patch(facecolor='green', alpha=0.15, label='Typical residual'),
                  Patch(facecolor='orange', alpha=0.15, label='Burst events')]
all_handles = line_elements + [intrinsic_line] + patch_elements
all_labels = [h.get_label() for h in all_handles]
ax.legend(all_handles, all_labels, loc='upper left', fontsize=11, framealpha=0.9)

ax.set_xlim(1e-10, 1e-3)
ax.set_ylim(1e-2, 1e5)
ax_T1.set_ylim(1e-2, 1e5)

plt.tight_layout()
plt.savefig(f'{dir_t1_calculation}/gamma_qp_dual_axis.png', dpi=300, bbox_inches='tight')
print(f"Dual-axis plot saved to {dir_t1_calculation}")

#==============================================================================
# TIME-DEPENDENT T1 CALCULATIONS (if data available)
#==============================================================================
if USE_TIME_DEPENDENT_DATA:
    print("\n" + "=" * 80)
    print("TIME-DEPENDENT T1 CALCULATIONS")
    print("=" * 80)

    Gamma_qp_t = prefactor * x_qp_t
    Gamma_total_t = Gamma_qp_t + Gamma_intrinsic

    T1_qp_t = 1 / Gamma_qp_t
    T1_total_t = 1 / Gamma_total_t

    Gamma_qp_t_kHz = Gamma_qp_t / 1e3
    Gamma_total_t_kHz = Gamma_total_t / 1e3
    T1_qp_t_us = T1_qp_t * 1e6
    T1_total_t_us = T1_total_t * 1e6

    t_us: np.ndarray = t_array * 1e6
    idx_peak_Gamma = np.argmax(Gamma_qp_t)
    idx_min_T1_qp = np.argmin(T1_qp_t)
    idx_min_T1_total = np.argmin(T1_total_t)

    print(f"\nTime-dependent results:")
    print(f"  Peak Gamma_qp: {Gamma_qp_t[idx_peak_Gamma]/1e3:.2f} kHz at t = {t_us[idx_peak_Gamma]:.2f} us")
    print(f"  Minimum T1_qp (QP-only): {T1_qp_t_us[idx_min_T1_qp]:.2f} us at t = {t_us[idx_min_T1_qp]:.2f} us")
    print(f"  Minimum T1_total (QP + intrinsic): {T1_total_t_us[idx_min_T1_total]:.2f} us at t = {t_us[idx_min_T1_total]:.2f} us")
    print(f"  Final Gamma_qp: {Gamma_qp_t[-1]/1e3:.2f} kHz")
    print(f"  Final T1_qp (QP-only): {T1_qp_t_us[-1]:.2f} us")
    print(f"  Final T1_total (QP + intrinsic): {T1_total_t_us[-1]:.2f} us")

    Gamma_qp_equilibrium = prefactor * x_qp_t[-1]
    T1_total_equilibrium = 1 / (Gamma_qp_equilibrium + Gamma_intrinsic)
    T1_total_equilibrium_us = T1_total_equilibrium * 1e6

    print(f"\n" + "=" * 80)
    print("T1 RECOVERY TIME ANALYSIS")
    print("=" * 80)
    print(f"Equilibrium T1_total: {T1_total_equilibrium_us:.2f} us")
    print(f"  (corresponding to final x_qp = {x_qp_t[-1]:.2e})")
    print(f"  (expected equilibrium x_qp_eq from constants = {x_qp_eq:.2e})")

    recovery_thresholds = [0.90, 0.95, 0.99]

    print(f"\nRecovery times to various thresholds:")
    for threshold in recovery_thresholds:
        target_T1 = threshold * T1_total_equilibrium_us

        if T1_total_t_us[idx_min_T1_total] >= target_T1:
            print(f"  {threshold*100:.0f}% recovery ({target_T1:.2f} us): T1 never dropped below this threshold")
            continue

        recovery_idx = np.where(T1_total_t_us[idx_min_T1_total:] >= target_T1)[0]

        if len(recovery_idx) > 0:
            recovery_time_idx = idx_min_T1_total + recovery_idx[0]
            recovery_time_us = t_us[recovery_time_idx]

            if recovery_time_us < 1e3:
                print(f"  {threshold*100:.0f}% recovery ({target_T1:.2f} us): {recovery_time_us:.2f} us")
            elif recovery_time_us < 1e6:
                print(f"  {threshold*100:.0f}% recovery ({target_T1:.2f} us): {recovery_time_us/1e3:.2f} ms")
            else:
                print(f"  {threshold*100:.0f}% recovery ({target_T1:.2f} us): {recovery_time_us/1e6:.2f} s")
        else:
            print(f"  {threshold*100:.0f}% recovery ({target_T1:.2f} us): Not yet recovered within simulation time")

    recovery_percentage = (T1_total_t_us[-1] - T1_total_t_us[idx_min_T1_total]) / (T1_total_equilibrium_us - T1_total_t_us[idx_min_T1_total]) * 100

    print(f"\nRecovery status at end of simulation:")
    print(f"  Initial T1 (at t=0): {T1_total_t_us[0]:.2f} us")
    print(f"  Minimum T1 (during burst): {T1_total_t_us[idx_min_T1_total]:.2f} us")
    print(f"  Final T1 (at t={t_array[-1]*1e3:.1f} ms): {T1_total_t_us[-1]:.2f} us")
    print(f"  Expected equilibrium T1: {T1_total_equilibrium_us:.2f} us")
    print(f"")
    print(f"  Recovery progress: {recovery_percentage:.1f}% of full recovery")
    if recovery_percentage >= 99.9:
        print(f"  STATUS: T1 has FULLY RECOVERED to equilibrium!")
    elif recovery_percentage >= 95:
        print(f"  STATUS: T1 has substantially recovered (>95%)")
        print(f"  Still {T1_total_equilibrium_us - T1_total_t_us[-1]:.2f} us away from equilibrium")
    else:
        print(f"  STATUS: T1 has NOT fully recovered")
        print(f"  Still {T1_total_equilibrium_us - T1_total_t_us[-1]:.2f} us away from equilibrium")
    
    print("=" * 80)

    # Create a combined plot showing correlation between x_qp(t) and T1(t)
    fig4, ax1 = plt.subplots(1, 1, figsize=(14, 8))
    ax2 = ax1.twinx()

    # Plot x_qp(t) on left axis
    line1, = ax1.semilogx(t_us, x_qp_t, 'b-', linewidth=3, label='x_qp(t)')
    ax1.axhline(x_qp_eq, color='blue', linestyle='--', linewidth=1.5, alpha=0.5)

    # Plot T1(t) on right axis
    line2a, = ax2.loglog(t_us, T1_qp_t_us, 'r--', linewidth=2, alpha=0.5, label='T1_qp(t) (QP-only)')
    line2b, = ax2.loglog(t_us, T1_total_t_us, 'r-', linewidth=3, label='T1_total(t)')
    ax2.axhline(T1_intrinsic * 1e6, color='purple', linestyle='--', linewidth=2, alpha=0.5)

    # Mark key points
    ax1.plot(t_us[idx_peak_Gamma], x_qp_t[idx_peak_Gamma], 'ko', markersize=12, zorder=5)
    ax2.plot(t_us[idx_min_T1_total], T1_total_t_us[idx_min_T1_total], 'ko', markersize=12, zorder=5)

    # Labels
    ax1.set_xlabel('Time after radiation event (us)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Normalized QP Density x_qp', fontsize=14, fontweight='bold', color='b')
    ax2.set_ylabel('Coherence Time T1 (us)', fontsize=14, fontweight='bold', color='r')

    ax1.set_title('Qubit Coherence Degradation from Cosmic Ray Impact\n' +
                  'RT Equations + Phonon Injection -> T1 Calculation',
                  fontsize=15, fontweight='bold', pad=20)

    # Color the tick labels
    ax1.tick_params(axis='y', labelcolor='b', labelsize=12)
    ax2.tick_params(axis='y', labelcolor='r', labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)

    # Grid
    ax1.grid(True, alpha=0.3, which='both')

    # Combined legend
    lines = [line1, line2a, line2b]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left', fontsize=13, framealpha=0.9)

    # Annotation
    ax1.annotate(f'Minimum T1_total = {T1_total_t_us[idx_min_T1_total]:.1f} us at t = {t_us[idx_min_T1_total]:.2f} us',
                xy=(t_us[idx_min_T1_total], x_qp_t[idx_peak_Gamma]),
                xytext=(t_us[idx_min_T1_total]*3, x_qp_t[idx_peak_Gamma]),
                fontsize=12, fontweight='bold',
                arrowprops=dict(arrowstyle='->', lw=2.5, color='black'),
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

    ax1.set_xlim(left=max(t_us[0], 1e-3))

    plt.tight_layout()
    plt.savefig(f'{dir_t1_calculation}/xqp_and_T1_vs_time.png', dpi=300, bbox_inches='tight')
    print(f"Combined x_qp(t) and T1(t) plot saved to: {dir_t1_calculation}/xqp_and_T1_vs_time.png")

    # Create plot with time in ms, gamma_qp in us^-1, and x_qp on right axis
    fig5, ax_gamma = plt.subplots(1, 1, figsize=(14, 8))
    ax_xqp = ax_gamma.twinx()

    t_ms: np.ndarray = t_array * 1e3  # Convert time to milliseconds
    Gamma_qp_t_us_inv = Gamma_qp_t / 1e6  # Convert Gamma_qp from s^-1 to us^-1

    # Plot Gamma_qp on left axis
    line1, = ax_gamma.semilogy(t_ms, Gamma_qp_t_us_inv, 'r-', linewidth=3, label='Gamma_qp(t)')
    ax_gamma.plot(t_ms[idx_peak_Gamma], Gamma_qp_t_us_inv[idx_peak_Gamma], 'ko',
                  markersize=12, zorder=5)

    # Plot x_qp on right axis
    line2, = ax_xqp.semilogy(t_ms, x_qp_t, 'b-', linewidth=3, label='x_qp(t)')
    ax_xqp.axhline(x_qp_eq, color='blue', linestyle='--', linewidth=1.5, alpha=0.5,
                   label=f'x_qp equilibrium')
    ax_xqp.plot(t_ms[idx_peak_Gamma], x_qp_t[idx_peak_Gamma], 'ko', markersize=12, zorder=5)

    # Labels and formatting
    ax_gamma.set_xlabel('Time (ms)', fontsize=14, fontweight='bold')
    ax_gamma.set_ylabel('QP Decay Rate Gamma_qp (μs⁻¹)', fontsize=14, fontweight='bold', color='r')
    ax_xqp.set_ylabel('Normalized QP Density x_qp', fontsize=14, fontweight='bold', color='b')

    ax_gamma.set_title('Quasiparticle Decay Rate and Density vs Time\n' +
                       'Following Cosmic Ray Impact on Transmon Qubit',
                       fontsize=15, fontweight='bold', pad=20)

    # Color the tick labels
    ax_gamma.tick_params(axis='y', labelcolor='r', labelsize=12)
    ax_xqp.tick_params(axis='y', labelcolor='b', labelsize=12)
    ax_gamma.tick_params(axis='x', labelsize=12)

    # Grid
    ax_gamma.grid(True, alpha=0.3, which='both')

    # Combined legend
    lines = [line1, line2]
    labels = [l.get_label() for l in lines]
    ax_gamma.legend(lines, labels, loc='upper right', fontsize=13, framealpha=0.9)

    # Annotation at peak
    ax_gamma.annotate(f'Peak: Gamma_qp = {Gamma_qp_t_us_inv[idx_peak_Gamma]:.3f} μs⁻¹\n' +
                      f'x_qp = {x_qp_t[idx_peak_Gamma]:.2e}\n' +
                      f't = {t_ms[idx_peak_Gamma]:.4f} ms',
                      xy=(t_ms[idx_peak_Gamma], Gamma_qp_t_us_inv[idx_peak_Gamma]),
                      xytext=(t_ms[idx_peak_Gamma] + 10, Gamma_qp_t_us_inv[idx_peak_Gamma] * 0.5),
                      fontsize=11, fontweight='bold',
                      arrowprops=dict(arrowstyle='->', lw=2, color='black'),
                      bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

    ax_gamma.set_xlim(0, t_ms[-1] * 0.25) # 25% of time range

    plt.tight_layout()
    plt.savefig(f'{dir_t1_calculation}/gamma_qp_and_xqp_vs_time_ms.png', dpi=300, bbox_inches='tight')
    print(f"Gamma_qp and x_qp vs time (ms) plot saved to: {dir_t1_calculation}/gamma_qp_and_xqp_vs_time_ms.png")

plt.show()

#==============================================================================
# PAIR-BREAKING EFFICIENCY CALCULATION
#==============================================================================
if USE_TIME_DEPENDENT_DATA:
    print("\n" + "=" * 80)
    print("PAIR-BREAKING EFFICIENCY (Energy in Al -> Energy in QPs)")
    print("=" * 80)

    # Load phonon injection data to get energy received by Al
    phonon_data_file = 'phonon_dynamics/phonon_injection_data.pkl'
    try:
        with open(phonon_data_file, 'rb') as f:
            phonon_data = pickle.load(f)

        E_received_Al = phonon_data['energy_received']  # Energy that reached Al film [J]

        # Calculate energy stored in excess quasiparticles
        # Use peak QP density for maximum excess QPs
        n_qp_peak = np.max(n_qp_t)
        n_qp_equilibrium = n_qp_t[-1]  # Use final value as equilibrium
        n_qp_excess = n_qp_peak - n_qp_equilibrium

        # Get Al volume from metadata
        V_Al = metadata['V_film']

        # Each QP carries energy ~ Delta
        Delta_Al_eV = metadata['Delta_Al']  # This is in eV
        Delta_Al_J = Delta_Al_eV * 1.602e-19  # Convert eV to Joules
        E_in_qps = n_qp_excess * V_Al * Delta_Al_J

        # Calculate pair-breaking efficiency
        eta_pair_breaking = (E_in_qps / E_received_Al) * 100

        print(f"\nEnergy received by Al film: {E_received_Al:.3e} J ({E_received_Al/1.602e-19:.2f} eV)")
        print(f"Peak QP density: {n_qp_peak:.3e} m^-3")
        print(f"Equilibrium QP density: {n_qp_equilibrium:.3e} m^-3")
        print(f"Excess QP density: {n_qp_excess:.3e} m^-3")
        print(f"Al film volume: {V_Al:.3e} m^3")
        print(f"Energy per QP (Delta): {Delta_Al_J:.3e} J ({Delta_Al_eV*1e6:.1f} ueV)")
        print(f"Energy stored in excess QPs: {E_in_qps:.3e} J ({E_in_qps/1.602e-19:.2f} eV)")
        print(f"\n*** PAIR-BREAKING EFFICIENCY: {eta_pair_breaking:.2f}% ***")
        print(f"    (Compare to Martinis 2021: ~57% for phonons with E > 2Delta)")
        print("=" * 80)

    except FileNotFoundError:
        print(f"\nCould not calculate pair-breaking efficiency: {phonon_data_file} not found")

print("\nAnalysis complete!")
