"""
Ballistic Phonon Transport in Superconducting Qubits

This script calculates ballistic phonon propagation from a point source
in Si substrate to a qubit location on the surface.

Key Physics:
- Phonons propagate ballistically (straight lines) from source
- Arrive at different times based on distance
- Geometric dilution follows 1/r³ scaling
- Interface transmission from Si to Al

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pickle
import os
import sys

# Add parent directory to path for importing constants
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from constants import (
    e, k_B, eV_to_J,
    Delta_Al, Delta_Al_eV, E_phonon,
    d_Al, A_qubit, V_Al,
    Si_sound_velocity as v_sound_Si, Z_Si, Z_Al, T_interface,
    x_0, y_0, z_0,
    chip_size, grid_resolution,
    t_max_phonon, n_points_phonon,
    time_spread_phonon,
    file_energy_deposit, dir_phonon_dynamics
)

print("="*80)
print("BALLISTIC PHONON TRANSPORT DEMONSTRATION")
print("="*80)

# ============================================================================
# LOAD ENERGY DEPOSIT DATA FROM BETHE-BLOCH CALCULATION
# ============================================================================
print("\n" + "="*80)
print("LOADING ENERGY DEPOSIT DATA")
print("="*80)

energy_deposit_file = file_energy_deposit
try:
    with open(energy_deposit_file, 'rb') as f:
        energy_data = pickle.load(f)

    E_deposited = energy_data['energy_deposited_joules']
    z_0 = energy_data['substrate_thickness_um'] * 1e-6

    print(f"Loaded energy deposit data from: {energy_deposit_file}")
    print(f"  Muon energy: {energy_data['muon_energy_mev']:.1f} MeV ({energy_data['muon_energy_mev']/1000:.2f} GeV)")
    print(f"  Substrate thickness: {energy_data['substrate_thickness_um']:.1f} um")
    print(f"  Incident angle: {energy_data['incident_angle_deg']:.1f} degrees")
    print(f"  Energy deposited in Si: {E_deposited:.6e} J ({E_deposited/e:.6e} eV)")
    print(f"  Source depth (z_0): {z_0*1e6:.1f} um")

    USE_BETHE_DATA = True

except FileNotFoundError:
    print(f"X Energy deposit file not found: {energy_deposit_file}")
    print("  Run 'Energy_Deposit/bethe.py' first to generate energy deposit data")
    print("  Using default values for demonstration...")
    E_deposited = 1e6 * e
    z_0 = 200e-6
    USE_BETHE_DATA = False

# ============================================================================
# PHYSICAL PARAMETERS (from centralized constants file)
# ============================================================================

# All physical parameters are now imported from constants.py
# Delta_Al, E_phonon, d_Al, A_qubit, V_Al
# v_sound_Si, Z_Si, Z_Al, T_interface

print(f"\nPhysical Parameters (from constants.py):")
print(f"  Al film thickness: {d_Al*1e9:.1f} nm")
print(f"  Qubit area: {A_qubit*1e12:.1f} um^2")
print(f"  Superconducting gap: {Delta_Al_eV*1e6:.1f} ueV")
print(f"  Si sound velocity: {v_sound_Si:.0f} m/s")
print(f"  Interface transmission: {T_interface*100:.2f}%")

# ============================================================================
# RADIATION EVENT PARAMETERS (from centralized constants file)
# ============================================================================


print(f"\n{'='*80}")
print("RADIATION EVENT")
print(f"{'='*80}")
print(f"Source location: ({x_0*1e3:.2f}, {y_0*1e3:.2f}) mm at depth {z_0*1e6:.1f} um")
print(f"Energy deposited: {E_deposited/e*1e-6:.2f} MeV ({E_deposited:.6e} J)")
if USE_BETHE_DATA:
    print(f"(From Bethe-Bloch calculation: {energy_data['muon_energy_mev']/1000:.2f} GeV muon)")

# Define single qubit position
qubit_positions = [
    (1.2e-3, 1.2e-3, "Qubit", "orange"),
]

# ============================================================================
# BALLISTIC PHONON CALCULATION
# ============================================================================

def calculate_ballistic_phonon_arrival(x_q, y_q, x_0, y_0, z_0, E_0, t_array):
    """
    Calculate time-dependent phonon arrival at qubit location.

    This calculates ballistic phonon transport from point source to qubit.

    Parameters:
    -----------
    x_q, y_q : float
        Qubit position on surface (m)
    x_0, y_0, z_0 : float
        Point source position (m)
    E_0 : float
        Deposited energy (J)
    t_array : ndarray
        Time points for evaluation (s)

    Returns:
    --------
    phonon_flux_t : ndarray
        Phonon flux [m⁻³ s⁻¹] as function of time
    t_arrival : float
        Peak arrival time (s)
    distance : float
        Source-to-qubit distance (m)
    """

    # Distance from source to qubit
    distance = np.sqrt((x_q - x_0)**2 + (y_q - y_0)**2 + z_0**2)

    # Ballistic arrival time
    t_arrival = distance / v_sound_Si

    # Solid angle factor for geometric dilution
    # cos(θ) = z_0 / r, differential solid angle dΩ ~ A_qubit / r²
    solid_angle_factor = z_0 / (4 * np.pi * distance**3)

    # Total number of phonons transmitted into this qubit
    # Energy fraction reaching qubit area from point source
    E_at_qubit = T_interface * E_0 * solid_angle_factor * A_qubit

    # Number of phonons (assuming characteristic energy E_phonon)
    N_phonons_total = E_at_qubit / E_phonon

    # Phonon density in Al film volume
    n_phonons_total = N_phonons_total / V_Al

    # Time spread of phonon pulse
    # For ballistic transport, pulse width is limited by angular spread
    # from finite source size and phonon frequency bandwidth
    time_spread = time_spread_phonon

    # Phonon flux: Gaussian pulse centered at t_arrival
    # Units: [m⁻³ s⁻¹] - phonons arriving per unit volume per unit time
    gaussian_pulse = np.exp(-(t_array - t_arrival)**2 / (2 * time_spread**2))
    gaussian_pulse /= (time_spread * np.sqrt(2 * np.pi))

    phonon_flux_t = n_phonons_total * gaussian_pulse

    flux_peak_analytical = n_phonons_total / (time_spread * np.sqrt(2 * np.pi))

    return phonon_flux_t, t_arrival, distance, E_at_qubit, n_phonons_total, flux_peak_analytical

# ============================================================================
# CALCULATE PHONON ARRIVAL FOR EACH QUBIT
# ============================================================================

print(f"\n{'='*80}")
print("BALLISTIC TRANSPORT CALCULATIONS")
print(f"{'='*80}")

# Time array for simulation with dense sampling around phonon arrival times
t_max = t_max_phonon
t_eval = np.linspace(0, t_max, n_points_phonon)

results = []

for x_q, y_q, label, color in qubit_positions:
    print(f"\n{label}:")
    print(f"  Position: ({x_q*1e3:.2f}, {y_q*1e3:.2f}) mm")

    # Calculate ballistic phonon arrival
    phonon_flux_t, t_arrival, distance, E_at_qubit, n_phonons_total, flux_peak_analytical = \
        calculate_ballistic_phonon_arrival(x_q, y_q, x_0, y_0, z_0,
                                             E_deposited, t_eval)

    print(f"  Distance from source: {distance*1e3:.2f} mm")
    print(f"  Arrival time: {t_arrival*1e9:.2f} ns")
    print(f"  Energy received: {E_at_qubit/e:.2e} eV = {E_at_qubit/e*1e-3:.2e} keV")
    print(f"  Transmission efficiency: {E_at_qubit/E_deposited*100:.4f}%")
    print(f"  Total phonon density: {n_phonons_total*1e-18:.2e} um^-3")
    print(f"  Peak flux (analytical): {flux_peak_analytical:.2e} m^-3 s^-1")

    # Find peak values
    peak_flux = np.max(phonon_flux_t)
    idx_peak = np.argmax(phonon_flux_t)
    t_peak = t_eval[idx_peak]

    print(f"  Peak flux (numerical): {peak_flux:.2e} m^-3 s^-1")
    print(f"  Time to peak: {t_peak*1e9:.2f} ns")

    results.append({
        'label': label,
        'color': color,
        'x': x_q,
        'y': y_q,
        'distance': distance,
        't_arrival': t_arrival,
        'E_received': E_at_qubit,
        'n_phonons': n_phonons_total,
        'phonon_flux_t': phonon_flux_t,
        'peak_flux': peak_flux,
        't_peak': t_peak,
    })

# ============================================================================
# ENERGY CONSERVATION CHECK
# ============================================================================

print(f"\n{'='*80}")
print("ENERGY CONSERVATION CHECK")
print(f"{'='*80}")

total_energy_to_qubits = sum([r['E_received'] for r in results])
print(f"Total energy deposited: {E_deposited/e*1e-3:.2f} keV")
print(f"Total energy to the qubit: {total_energy_to_qubits/e*1e-3:.3f} keV")
print(f"Fraction captured: {total_energy_to_qubits/E_deposited*100:.4f}%")
print(f"Energy lost to substrate: {(1 - total_energy_to_qubits/E_deposited)*100:.2f}%")

# ============================================================================
# HEATMAP DATA CALCULATION
# ============================================================================

# Create a 2D grid for the heatmap
grid_res = grid_resolution
x_grid = np.linspace(0, chip_size, grid_res)
y_grid = np.linspace(0, chip_size, grid_res)
X, Y = np.meshgrid(x_grid, y_grid)

# Calculate phonon density across the grid
distance_grid = np.sqrt((X - x_0)**2 + (Y - y_0)**2 + z_0**2)
solid_angle_grid = z_0 / (4 * np.pi * distance_grid**3)
E_at_point_grid = T_interface * E_deposited * solid_angle_grid * A_qubit
N_phonons_grid = E_at_point_grid / E_phonon
phonon_density_grid = N_phonons_grid / V_Al

phonon_density_grid_um = phonon_density_grid * 1e-18

# ============================================================================
# VISUALIZATION
# ============================================================================

print(f"\n{'='*80}")
print("GENERATING PLOTS")
print(f"{'='*80}")

# Convert time to nanoseconds for plotting
t_ns = t_eval * 1e9

# Figure: Create a composite layout with time series and heatmap
fig = plt.figure(figsize=(18, 9))
gs = fig.add_gridspec(2, 2, width_ratios=(1, 1.2))

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[:, 1])

# Top-left: Phonon flux
for res in results:
    flux_per_ns = res['phonon_flux_t'] * 1e-9
    ax1.plot(t_ns, flux_per_ns, linewidth=2.5,
             label=f"{res['label']} (d={res['distance']*1e3:.1f} mm)",
             color=res['color'])
    ax1.axvline(res['t_arrival']*1e9, color=res['color'],
                linestyle=':', alpha=0.4, linewidth=1.5)

ax1.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Phonon Flux (m^-3 ns^-1)', fontsize=12, fontweight='bold')
ax1.set_title('Ballistic Phonon Arrival at the Qubit',
              fontsize=13, fontweight='bold')
ax1.legend(fontsize=10, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')
ax1.set_ylim([1e6, 1e21])

# Bottom-left: Cumulative energy arrived
for res in results:
    cumulative_energy = np.cumsum(res['phonon_flux_t']) * (t_eval[1] - t_eval[0]) * V_Al * E_phonon
    ax2.plot(t_ns, cumulative_energy/e*1e-3, linewidth=2.5,
             label=res['label'], color=res['color'])
    ax2.axhline(res['E_received']/e*1e-3, color=res['color'],
                linestyle='--', alpha=0.4, linewidth=1.5)

ax2.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Cumulative Energy (keV)', fontsize=12, fontweight='bold')
ax2.set_title('Energy Conservation: Integrated Phonon Flux',
              fontsize=13, fontweight='bold')
ax2.legend(fontsize=10, loc='lower right')
ax2.grid(True, alpha=0.3)

# Right: 2D Heatmap of Phonon Density
im = ax3.pcolormesh(X * 1e3, Y * 1e3, np.log10(phonon_density_grid_um + 1e-9),
                    shading='auto', cmap='inferno')
ax3.set_aspect('equal')

# Mark source and qubit locations
ax3.scatter(x_0 * 1e3, y_0 * 1e3, marker='*', s=300, edgecolor='white',
            facecolor='cyan', linewidth=1.5, label='Muon Hit (Source)')
for res in results:
    ax3.scatter(res['x'] * 1e3, res['y'] * 1e3, marker='s', s=150, edgecolor='white',
                facecolor=res['color'], linewidth=1.5, label=res['label'])

ax3.set_xlabel('X position (mm)', fontsize=12, fontweight='bold')
ax3.set_ylabel('Y position (mm)', fontsize=12, fontweight='bold')
ax3.set_title('Spatial Distribution of Peak Phonon Density',
              fontsize=13, fontweight='bold')
ax3.legend(fontsize=10)
fig.colorbar(im, ax=ax3, label='Log10(Peak Phonon Density [um^-3])')


plt.tight_layout(pad=3.0)
plt.savefig('Phonon_Dynamics/ballistic_phonon_flux_and_heatmap.png', dpi=300, bbox_inches='tight')
print("Saved: ballistic_phonon_flux_and_heatmap.png")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print(f"\n{'='*80}")
print("SUMMARY: BALLISTIC PHONON TRANSPORT")
print(f"{'='*80}")

# Get data for summary from results
arrivals = [r['t_arrival']*1e9 for r in results]
distances = [r['distance']*1e3 for r in results]
peak_fluxes_per_ns = [r['peak_flux']*1e-9 for r in results]

print(f"\nKey Timescales:")
print(f"  Phonon transit time: {arrivals[0]:.1f} ns")

print(f"\nQubit Properties:")
print(f"  Distance: {distances[0]:.2f} mm -> Peak flux = {peak_fluxes_per_ns[0]:.2e} m^-3 ns^-1")

print(f"\n{'='*80}")
print("CREATING PHONON INJECTION FUNCTION FOR RT EQUATIONS")
print(f"{'='*80}")

# ============================================================================
# CREATE PHONON INJECTION FUNCTION FOR RT EQUATIONS
# ============================================================================

# Get the phonon flux data
phonon_flux_data = results[0]['phonon_flux_t']

# Create interpolation function for phonon injection rate
phonon_injection_func = interp1d(
    t_eval,
    phonon_flux_data,
    kind='cubic',
    bounds_error=False,
    fill_value=0.0
)

# Save the function and metadata to a pickle file
phonon_data = {
    'phonon_injection_func': phonon_injection_func,
    't_array': t_eval,
    'phonon_flux': phonon_flux_data,
    't_arrival': results[0]['t_arrival'],
    'peak_flux': results[0]['peak_flux'],
    'total_phonon_density': results[0]['n_phonons'],
    'energy_received': results[0]['E_received'],
    'metadata': {
        'source_position': (x_0, y_0, z_0),
        'qubit_position': (results[0]['x'], results[0]['y'], 0),
        'distance': results[0]['distance'],
        'E_deposited': E_deposited,
        'E_phonon': E_phonon,
        'T_interface': T_interface,
        'v_sound_Si': v_sound_Si,
        'Al_film_thickness': d_Al,
        'qubit_area': A_qubit,
        'Al_volume': V_Al,
        'use_bethe_data': USE_BETHE_DATA,
        'bethe_muon_energy_mev': energy_data['muon_energy_mev'] if USE_BETHE_DATA else None,
    }
}

os.makedirs(dir_phonon_dynamics, exist_ok=True)
output_file = f'{dir_phonon_dynamics}/phonon_injection_data.pkl'
with open(output_file, 'wb') as f:
    pickle.dump(phonon_data, f)

print(f"\nPhonon injection function saved to: {output_file}")
print(f"\nFunction properties:")
print(f"  Time range: {t_eval[0]*1e9:.2f} ns to {t_eval[-1]*1e6:.2f} us")
print(f"  Peak arrival: {results[0]['t_arrival']*1e9:.2f} ns")
print(f"  Peak flux: {results[0]['peak_flux']:.2e} m^-3 s^-1")
print(f"  Total phonon density: {results[0]['n_phonons']:.2e} m^-3")
print(f"\nThis function can be used in RT equations as:")
print(f"  dn_ph/dt += phonon_injection_func(t)")

print(f"\n{'='*80}")
print("SIMULATION COMPLETE")
print(f"{'='*80}")

plt.show()
