"""
Centralized Constants for Transmon Qubit Modeling

This file contains all physical constants, material properties, and experimental
parameters used across the four simulation modules:
- energy_deposit: Bethe-Bloch muon energy deposition
- phonon_dynamics: Ballistic phonon transport
- rt_equations: Rothwarf-Taylor equations
- t1_calculation: T1 coherence time calculations

All values are in SI units unless otherwise specified.
"""

import numpy as np
from scipy import constants

#==============================================================================
# FUNDAMENTAL PHYSICAL CONSTANTS
#==============================================================================
h = constants.h  # Planck constant [J*s]
hbar = constants.hbar  # Reduced Planck constant [J*s]
k_B = constants.k  # Boltzmann constant [J/K]
k_B_eV_K = 8.617e-5  # Boltzmann constant [eV/K]
c = constants.c  # Speed of light [m/s]
e = constants.e  # Elementary charge [C]
N_A = constants.N_A  # Avogadro's number [mol^-1]
m_e = constants.m_e  # Electron mass [kg]
m_e_MeV = 0.510998946  # Electron mass [MeV/c^2]

#==============================================================================
# UNIT CONVERSIONS
#==============================================================================
eV_to_J = constants.e  # Conversion factor from eV to Joules
J_to_eV = 1 / eV_to_J  # Conversion factor from Joules to eV
MeV_to_J = 1.602176634e-13  # Conversion factor from MeV to Joules

#==============================================================================
# ALUMINUM SUPERCONDUCTOR PARAMETERS
#==============================================================================
# Superconducting gap
Delta_Al_eV = 180e-6  # Superconducting gap [eV]
Delta_Al = Delta_Al_eV * eV_to_J  # Superconducting gap [J]

# Density of states at Fermi level (single-spin)
N_0_per_eV = 1.74e28  # Single-spin DOS at Fermi level [eV^-1 m^-3]
N_cp = 2 * N_0_per_eV * Delta_Al_eV  # Cooper pair density [m^-3]

# Material properties
Al_density = 2700  # Aluminum density [kg/m^3]
Al_sound_velocity = 3100  # Sound velocity in Al [m/s]
Z_Al = Al_density * Al_sound_velocity  # Acoustic impedance of Al [kg/(m^2*s)]

#==============================================================================
# ALUMINUM FILM GEOMETRY
#==============================================================================
d_Al = 50e-9  # Al film thickness [m] (50 nm)
A_qubit = (100e-6)**2  # Qubit area [m^2] (100 um x 100 um)
V_Al = d_Al * A_qubit  # Al film volume [m^3]

#==============================================================================
# SILICON SUBSTRATE PARAMETERS
#==============================================================================
# Material properties for phonon transport
Si_density = 2330  # Silicon density [kg/m^3]
Si_sound_velocity = 5500  # Sound velocity in Si [m/s]
Z_Si = Si_density * Si_sound_velocity  # Acoustic impedance of Si [kg/(m^2*s)]

# Material properties for Bethe-Bloch calculations
Si_Z = 14  # Atomic number
Si_A = 28.0855  # Atomic mass [g/mol]
Si_density_g_cm3 = 2.329  # Density [g/cm^3]
Si_I_eV = 173.0  # Mean excitation energy [eV]

# Sternheimer density correction parameters for Silicon
Si_x0 = 0.2014
Si_x1 = 2.8716
Si_C = -4.4355
Si_a = 0.14921
Si_m = 3.2546

# Interface transmission coefficient (Si → Al)
T_interface = 4 * Z_Si * Z_Al / (Z_Si + Z_Al)**2

#==============================================================================
# TRANSMON QUBIT PARAMETERS
#==============================================================================
f_q = 6e9  # Qubit frequency [Hz] (6 GHz typical)
omega_q = 2 * np.pi * f_q  # Angular frequency [rad/s]

#==============================================================================
# MUON/COSMIC RAY PARAMETERS
#==============================================================================
muon_mass_MeV = 105.6583745  # Muon rest mass [MeV/c^2]
muon_charge = 1  # Muon charge in units of e
muon_energy_MeV = 2000  # Typical muon kinetic energy [MeV] (2 GeV)
substrate_thickness_um = 200  # Silicon substrate thickness [um]
incident_angle_deg = 0  # Muon incident angle from normal [degrees]

# Bethe-Bloch constant
K_BetheBloch = 0.307075  # [MeV*cm^2/g] for A=1 g/mol

#==============================================================================
# PHONON PARAMETERS
#==============================================================================
E_phonon = 3 * Delta_Al  # Characteristic phonon energy for pair breaking [J]
time_spread_phonon = 50e-9  # Phonon pulse width [s] (50 ns)

# Phonon timescales
tau_esc = 80e-9  # Phonon escape time from Al film [s] (80 ns)
tau_anh = 1e-9  # Anharmonic decay time [s] (1 ns)

#==============================================================================
# ROTHWARF-TAYLOR EQUATION PARAMETERS
#==============================================================================
# Temperature
T_bath = 0.020  # Bath temperature [K] (20 mK)

# Rate coefficients
R = 1.0e-18  # Recombination coefficient [m^3/s]
B = 9e9  # Pair-breaking rate [s^-1] (1/tau_pb from literature for Al)

# Quasiparticle parameters
tau_trap = 1e-3  # Quasiparticle trapping/relaxation time [s] (1 ms)
                  # Based on experimental measurements: Nature Comm 4, 5836 (2014)
                  # Reports ~0.794 ms tunneling time and ms-scale T1 recovery
f_pb = 0.5  # Fraction of recombination phonons with E > 2Δ (empirical)

# Equilibrium quasiparticle density (experimental values)
# Measured: 25-55 per μm³ below 160 mK
n_qp_eq = 1e18  # Equilibrium QP density [m^-3] (1 per μm³)
x_qp_eq = n_qp_eq / (2 * N_0_per_eV * Delta_Al_eV)  # Normalized equilibrium x_qp

# Equilibrium phonon density (E > 2Δ)
n_ph_eq = 0.0  # Equilibrium high-energy phonon density [m^-3]

#==============================================================================
# SIMULATION PARAMETERS
#==============================================================================
# Energy deposition
use_bethe_data = True  # Use Bethe-Bloch calculation for energy deposit

# Phonon transport
chip_size = 4e-3  # Chip size for spatial grid [m] (4 mm)
grid_resolution = 200  # Grid resolution for heatmap plots

# Source location on chip (for phonon transport)
x_0 = 0.5e-3  # Source x position [m] (0.5 mm)
y_0 = 0.5e-3  # Source y position [m] (0.5 mm)
z_0 = 200e-6  # Source depth in substrate [m] (200 um) - updated from Bethe-Bloch

# Time arrays
t_max_phonon = 2e-6  # Max time for phonon simulation [s] (2 us)
n_points_phonon = 5000  # Number of time points for phonon simulation

t_max_RT = 100e-3  # Max time for RT equations [s] (1000 ms)
n_points_RT = 2000  # Number of time points for RT equations

#==============================================================================
# T1 CALCULATION PARAMETERS
#==============================================================================
# Prefactor for Gamma_qp calculation
# Gamma_qp = sqrt(2*omega_q*Delta/(pi^2*hbar)) * x_qp
Gamma_qp_prefactor = np.sqrt(2 * omega_q * Delta_Al / (np.pi**2 * hbar))  # [s^-1]

# Intrinsic T1 (non-QP loss mechanisms)
T1_intrinsic = 100e-6  # Intrinsic coherence time [s] (100 μs)
Gamma_intrinsic = 1 / T1_intrinsic  # Intrinsic decay rate [s^-1]

# Typical x_qp values for reference
x_qp_typical_residual = 1e-7  # Typical residual background
x_qp_burst_low = 1e-6  # Low-level burst
x_qp_burst_high = 1e-5  # High-level burst

#==============================================================================
# OUTPUT DIRECTORIES
#==============================================================================
dir_energy_deposit = 'energy_deposit'
dir_phonon_dynamics = 'phonon_dynamics'
dir_rt_equations = 'rt_equations'
dir_t1_calculation = 't1_calculation'

#==============================================================================
# DATA FILE NAMES
#==============================================================================
file_energy_deposit = 'energy_deposit/energy_deposit_data.pkl'
file_phonon_injection = 'phonon_dynamics/phonon_injection_data.pkl'
file_xqp_vs_time = 'rt_equations/xqp_vs_time_data.pkl'

#==============================================================================
# PLOT CONFIGURATION
#==============================================================================
plot_dpi = 300  # DPI for saved plots
plot_figsize_single = (10, 8)  # Figure size for single plots
plot_figsize_multi = (14, 10)  # Figure size for multi-panel plots
