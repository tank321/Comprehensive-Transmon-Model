"""
Bethe-Bloch Energy Deposition Calculator for Muons in Silicon

Calculates energy deposited by cosmic ray muons in Si substrate and exports
the result for use in phonon transport simulations.

Physics: Muons traversing matter lose energy via ionization described by the
Bethe-Bloch formula. This energy is converted to phonons in the substrate.
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple
import matplotlib.pyplot as plt
import pickle
import os
import sys

# Add parent directory to path for importing constants
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from constants import (
    N_A, c, m_e_MeV, MeV_to_J, K_BetheBloch,
    muon_mass_MeV, muon_charge, muon_energy_MeV,
    substrate_thickness_um, incident_angle_deg,
    Si_Z, Si_A, Si_density_g_cm3, Si_I_eV,
    Si_x0, Si_x1, Si_C, Si_a, Si_m,
    dir_energy_deposit
)

@dataclass
class MaterialProperties:
    """Material properties for energy deposition calculations"""
    Z: float  # Atomic number
    A: float  # Atomic mass (g/mol)
    density: float  # g/cm^3
    I: float  # Mean excitation energy (eV)
    thickness: float  # Material thickness (cm)
    # Sternheimer density correction parameters
    x0: float = 0.2014
    x1: float = 2.8716
    C: float = -4.4355
    a: float = 0.14921
    m: float = 3.2546

@dataclass
class ParticleProperties:
    """Incident particle properties"""
    mass: float  # Rest mass (MeV/c^2)
    charge: int  # Charge in units of e
    energy: float  # Initial kinetic energy (MeV)

class BetheBlochSimulator:
    """
    Bethe-Bloch energy deposition calculator for muons in materials
    """

    def __init__(self):
        # Physical constants (from centralized constants file)
        self.m_e = m_e_MeV  # MeV/c^2 (electron mass)
        self.c = c * 1e2  # cm/s
        self.N_A = N_A  # mol^-1
        self.K = K_BetheBloch  # MeV*cm^2/g for A=1 g/mol

        # Muon properties (from centralized constants file)
        self.muon = ParticleProperties(
            mass=muon_mass_MeV,  # MeV/c^2
            charge=muon_charge,
            energy=0
        )

    def relativistic_params(self, kinetic_energy: float) -> Tuple[float, float]:
        """Calculate beta and gamma from kinetic energy"""
        gamma = 1 + kinetic_energy / self.muon.mass
        beta = np.sqrt(1 - 1/gamma**2) if gamma > 1 else 0
        return beta, gamma

    def max_energy_transfer(self, kinetic_energy: float) -> float:
        """Calculate maximum energy transfer T_max"""
        beta, gamma = self.relativistic_params(kinetic_energy)

        M = self.muon.mass
        m = self.m_e

        numerator = 2 * m * beta**2 * gamma**2
        denominator = 1 + 2*gamma*(m/M) + (m/M)**2

        T_max = numerator / denominator
        return min(T_max, kinetic_energy)

    def density_correction(self, beta_gamma: float, material: MaterialProperties) -> float:
        """Sternheimer density correction delta(beta*gamma)"""
        x = np.log10(beta_gamma)

        if x < material.x0:
            delta = 0
        elif x < material.x1:
            delta = 4.6052 * x + material.C + material.a * (material.x1 - x)**material.m
        else:
            delta = 4.6052 * x + material.C

        return max(delta, 0)

    def stopping_power(self, kinetic_energy: float, material: MaterialProperties) -> float:
        """
        Calculate mean stopping power -dE/dx using Bethe-Bloch formula
        Returns: MeV*cm^2/g
        """
        if kinetic_energy <= 0:
            return 0

        beta, gamma = self.relativistic_params(kinetic_energy)

        if beta < 1e-6:
            return 0

        T_max = self.max_energy_transfer(kinetic_energy)
        beta_gamma = beta * gamma
        delta = self.density_correction(beta_gamma, material)

        I_MeV = material.I * 1e-6

        # Bethe-Bloch formula
        prefactor = self.K * self.muon.charge**2 * (material.Z / material.A) / beta**2
        log_term = np.log(2 * self.m_e * beta**2 * gamma**2 * T_max / I_MeV**2)
        corrections = 2 * beta**2 + delta

        dE_dx = prefactor * (0.5 * log_term - corrections)

        return max(dE_dx, 0)

    def energy_deposit(self, kinetic_energy: float,
                      material: MaterialProperties,
                      angle_deg: float = 0,
                      use_integration: bool = True) -> float:
        """
        Calculate energy deposited in material

        Args:
            kinetic_energy: Initial muon kinetic energy (MeV)
            material: Material properties
            angle_deg: Incident angle from normal (degrees)
            use_integration: Use numerical integration for accurate calculation

        Returns:
            Energy deposited (MeV)
        """
        # Account for incident angle
        path_length = material.thickness / np.cos(np.radians(angle_deg))

        if use_integration or path_length > 0.01:
            # Use numerical integration for accuracy
            n_steps = max(100, int(path_length * 10000))
            dx = path_length / n_steps
            E_current = kinetic_energy
            E_deposited = 0

            for _ in range(n_steps):
                if E_current <= 0:
                    break
                dE = self.stopping_power(E_current, material) * material.density * dx
                dE = min(dE, E_current)
                E_deposited += dE
                E_current -= dE

            return E_deposited
        else:
            # Simple approximation for thin targets
            stopping_power = self.stopping_power(kinetic_energy, material)
            dE_dx_linear = stopping_power * material.density
            energy_deposited = dE_dx_linear * path_length
            return min(energy_deposited, kinetic_energy)

def calculate_energy_deposit(muon_energy_mev=muon_energy_MeV,
                            substrate_thickness_um=substrate_thickness_um,
                            incident_angle_deg=incident_angle_deg):
    """
    Main function to calculate energy deposited by muon in Si substrate

    Args:
        muon_energy_mev: Muon kinetic energy (MeV)
        substrate_thickness_um: Silicon substrate thickness (micrometers)
        incident_angle_deg: Incident angle from normal (degrees)

    Returns:
        dict: Energy deposit results
    """
    print("=" * 80)
    print("BETHE-BLOCH ENERGY DEPOSITION CALCULATION")
    print("=" * 80)

    # Create simulator
    sim = BetheBlochSimulator()

    # Silicon substrate properties (from centralized constants file)
    silicon = MaterialProperties(
        Z=Si_Z,
        A=Si_A,
        density=Si_density_g_cm3,  # g/cm^3
        I=Si_I_eV,  # eV
        thickness=substrate_thickness_um * 1e-4,  # Convert um to cm
        x0=Si_x0,
        x1=Si_x1,
        C=Si_C,
        a=Si_a,
        m=Si_m
    )

    print(f"\nSimulation Parameters:")
    print(f"  Muon kinetic energy: {muon_energy_mev} MeV ({muon_energy_mev/1000:.2f} GeV)")
    print(f"  Silicon thickness: {substrate_thickness_um} um")
    print(f"  Incident angle: {incident_angle_deg} degrees")
    print(f"  Path length: {silicon.thickness / np.cos(np.radians(incident_angle_deg)) * 1e4:.2f} um")

    # Calculate stopping power
    beta, gamma = sim.relativistic_params(muon_energy_mev)
    stopping_power = sim.stopping_power(muon_energy_mev, silicon)

    print(f"\nMuon Properties:")
    print(f"  Beta: {beta:.6f}")
    print(f"  Gamma: {gamma:.4f}")
    print(f"  Beta*Gamma: {beta*gamma:.4f}")

    print(f"\nStopping Power:")
    print(f"  dE/dx (mass): {stopping_power:.4f} MeV*cm^2/g")
    print(f"  dE/dx (linear): {stopping_power * silicon.density:.4f} MeV/cm")

    # Calculate energy deposit
    energy_deposited_mev = sim.energy_deposit(
        muon_energy_mev,
        silicon,
        incident_angle_deg,
        use_integration=True
    )

    # Convert to Joules (for phonon calculations) using centralized constant
    energy_deposited_joules = energy_deposited_mev * MeV_to_J

    print(f"\nEnergy Deposition Results:")
    print(f"  Energy deposited: {energy_deposited_mev:.6f} MeV")
    print(f"  Energy deposited: {energy_deposited_joules:.6e} J")
    print(f"  Fraction of muon energy: {energy_deposited_mev/muon_energy_mev*100:.4f}%")

    # Package results
    stochastic_samples_mev = None
    results = {
        'muon_energy_mev': muon_energy_mev,
        'substrate_thickness_um': substrate_thickness_um,
        'incident_angle_deg': incident_angle_deg,
        'path_length_cm': silicon.thickness / np.cos(np.radians(incident_angle_deg)),
        'stopping_power_mev_cm2_g': stopping_power,
        'energy_deposited_mev': energy_deposited_mev,
        'energy_deposited_joules': energy_deposited_joules,
        'beta': beta,
        'gamma': gamma,
        'stochastic_samples_mev': stochastic_samples_mev,
        'material_properties': {
            'Z': silicon.Z,
            'A': silicon.A,
            'density_g_cm3': silicon.density,
            'I_eV': silicon.I,
            'thickness_cm': silicon.thickness,
        }
    }

    return results

def save_energy_deposit(results, output_dir=dir_energy_deposit):
    """Save energy deposit results to file for use by phonon script"""

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'energy_deposit_data.pkl')

    with open(output_file, 'wb') as f:
        pickle.dump(results, f)

    print(f"\n{'='*80}")
    print("ENERGY DEPOSIT DATA SAVED")
    print(f"{'='*80}")
    print(f"Saved to: {output_file}")
    print(f"\nThis data contains:")
    print(f"  - Energy deposited in Si: {results['energy_deposited_joules']:.6e} J")
    print(f"  - Muon parameters (energy, angle)")
    print(f"  - Material properties")
    print(f"\nThis can be used as input for phonon transport simulations.")

    return output_file

def plot_energy_deposit_summary(results):
    """Create summary plots of energy deposition"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Create simulator for additional calculations
    sim = BetheBlochSimulator()
    silicon = MaterialProperties(
        Z=Si_Z, A=Si_A, density=Si_density_g_cm3, I=Si_I_eV,
        thickness=results['substrate_thickness_um'] * 1e-4,
        x0=Si_x0, x1=Si_x1, C=Si_C, a=Si_a, m=Si_m
    )

    # Plot 1: Stopping power vs energy
    ax = axes[0, 0]
    energies = np.logspace(2, 5, 200)  # 100 MeV to 100 GeV
    stopping_powers = [sim.stopping_power(E, silicon) for E in energies]
    ax.loglog(energies/1000, stopping_powers, 'b-', linewidth=2)
    ax.axvline(results['muon_energy_mev']/1000, color='r', linestyle='--',
               linewidth=2, label=f"This event: {results['muon_energy_mev']/1000:.1f} GeV")
    ax.axhline(1.664, color='green', linestyle='--', alpha=0.5, label='MIP value')
    ax.set_xlabel('Muon Energy (GeV)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Stopping Power (MeV*cm^2/g)', fontsize=12, fontweight='bold')
    ax.set_title('Bethe-Bloch Stopping Power in Silicon', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')

    # Plot 2: Energy deposit vs thickness
    ax = axes[0, 1]
    thicknesses_um = np.linspace(50, 500, 100)
    deposits = []
    for thick in thicknesses_um:
        silicon_temp = MaterialProperties(
            Z=Si_Z, A=Si_A, density=Si_density_g_cm3, I=Si_I_eV,
            thickness=thick * 1e-4,
            x0=Si_x0, x1=Si_x1, C=Si_C, a=Si_a, m=Si_m
        )
        dep = sim.energy_deposit(results['muon_energy_mev'], silicon_temp)
        deposits.append(dep)

    ax.plot(thicknesses_um, deposits, 'b-', linewidth=2)
    ax.axvline(results['substrate_thickness_um'], color='r', linestyle='--',
               linewidth=2, label=f"This substrate: {results['substrate_thickness_um']} um")
    ax.axhline(results['energy_deposited_mev'], color='r', linestyle='--',
               linewidth=1, alpha=0.5)
    ax.set_xlabel('Si Substrate Thickness (um)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Energy Deposited (MeV)', fontsize=12, fontweight='bold')
    ax.set_title(f'Energy Deposit vs Thickness ({results["muon_energy_mev"]/1000:.1f} GeV Muon)',
                 fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Energy deposit vs incident angle
    ax = axes[1, 0]
    angles = np.linspace(0, 75, 100)
    deposits_angle = []
    for angle in angles:
        dep = sim.energy_deposit(results['muon_energy_mev'], silicon, angle)
        deposits_angle.append(dep)

    ax.plot(angles, deposits_angle, 'b-', linewidth=2)
    ax.axvline(results['incident_angle_deg'], color='r', linestyle='--',
               linewidth=2, label=f"This event: {results['incident_angle_deg']} deg")
    ax.set_xlabel('Incident Angle (degrees)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Energy Deposited (MeV)', fontsize=12, fontweight='bold')
    ax.set_title('Energy Deposit vs Incident Angle', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Stochastic distribution (if available)
    ax = axes[1, 1]
    if results['stochastic_samples_mev'] is not None:
        samples = results['stochastic_samples_mev']
        ax.hist(samples, bins=50, alpha=0.7, density=True, edgecolor='black', linewidth=0.5)
        ax.axvline(np.mean(samples), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {np.mean(samples):.4f} MeV')
        ax.axvline(np.median(samples), color='green', linestyle='--', linewidth=2,
                   label=f'Median: {np.median(samples):.4f} MeV')
        ax.set_xlabel('Energy Deposited (MeV)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Probability Density', fontsize=12, fontweight='bold')
        ax.set_title('Stochastic Energy Deposit Distribution (Landau Fluctuations)',
                     fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'Deterministic Calculation Only\\nNo Stochastic Samples',
                ha='center', va='center', fontsize=14, transform=ax.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

    plt.tight_layout()
    plt.savefig('energy_deposit/energy_deposit_summary.png', dpi=300, bbox_inches='tight')
    print(f"\nSummary plot saved to: energy_deposit/energy_deposit_summary.png")

    return fig

if __name__ == "__main__":
    # Calculate energy deposition for typical cosmic ray muon
    # Parameters now come from centralized constants file
    results = calculate_energy_deposit()

    # Save results for phonon script
    save_energy_deposit(results)

    # Create summary plots
    plot_energy_deposit_summary(results)

    plt.show()

    print(f"\n{'='*80}")
    print("SIMULATION COMPLETE")
    print(f"{'='*80}")
