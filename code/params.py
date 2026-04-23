"""
weyl_hydro/params.py

Parameters and characteristic scales for the Weyl semimetal hydrodynamic
thermoelectric calculation.

All quantities are in natural units where hbar = k_B = 1, unless otherwise noted.
Energies in meV, temperatures in K (converted where needed), lengths in nm.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

# Fundamental constants in "natural" and SI units where convenient
HBAR_meV_ps  = 0.6582119569      # hbar in meV*ps
KB_meV_per_K = 0.08617333262     # k_B in meV/K
E_CHARGE     = 1.602176634e-19   # C
HBAR_SI      = 1.054571817e-34   # J*s
KB_SI        = 1.380649e-23      # J/K


@dataclass
class WeylParams:
    """
    Physical parameters for the two-node Weyl hydrodynamic model.

    Energy scales are in meV, temperatures in K, times in ps, lengths in nm,
    velocities in nm/ps.

    The default values are chosen to roughly correspond to a clean Weyl
    semimetal similar in spirit to WP2 / NbP in the hydrodynamic temperature
    window (T ~ 5-20 K).
    """

    # Fermi velocity (nm / ps).  v_F = 3e5 m/s -> 300 nm/ps
    vF: float = 300.0

    # Chemical potential (meV), measured from Weyl node
    mu: float = 5.0

    # Temperature (K)
    T: float = 10.0

    # Momentum-relaxing (elastic / impurity) time tau_mr (ps).
    # For WP2, the low-T tau is long; 10 ps is representative.
    tau_mr: float = 10.0

    # Inelastic (electron-electron) scattering time tau_ee (ps).
    # In hydrodynamic regime we require tau_ee < tau_mr.
    tau_ee: float = 1.0

    # Intervalley / chirality-flipping scattering time tau_5 (ps).
    # Typically tau_5 >> tau_ee, often >> tau_mr as well.
    tau_5: float = 50.0

    # Anomalous Hall coefficient controlled by chiral shift b (1/nm).
    # sigma_xy^AHE = e^2 b / (2 pi^2 hbar) in SI units.
    # Here we track it as a dimensionless knob in (0, 1) scaled to sigma0.
    b_shift: float = 0.1  # in 1/nm, along z.

    # Magnetic field (T).
    B: float = 0.0

    # Width of slab (nm).
    W: float = 500.0

    # Geometry-related
    no_slip: bool = True        # no-slip BC on u_x(+-W/2)
    absorb_chiral: bool = True  # n_5(+-W/2) = 0

    # Optional tilt (dimensionless) for future use
    tilt: float = 0.0


def thermodynamic_quantities(p: WeylParams):
    """
    Thermodynamic quantities for a two-node Weyl gas with chemical potential mu
    (same at both nodes) and temperature T.  Following standard relativistic-
    like expressions for a single 3D Weyl cone (the gas is relativistic with
    linear dispersion E = v_F |p|).

    We compute, per Weyl node (sum factor of 2 applied where relevant):
      n  : charge density
      s  : entropy density
      w  : enthalpy density w = eps + P
      chi_5 : chiral susceptibility (dn_5 / d mu_5 at mu_5 = 0)

    These are evaluated using standard Sommerfeld / full integrals for a single
    3D massless fermion gas in the grand canonical ensemble:
        n    = (1/(6 pi^2 (hbar vF)^3)) * (mu^3 + pi^2 mu T^2)       (low T)
        eps  = (1/(8 pi^2 (hbar vF)^3)) * (mu^4 + 2 pi^2 mu^2 T^2 + 7 pi^4 T^4 /15)
        P    = eps / 3 (relativistic)

    For the purposes of the linearized hydrodynamic calculation we need
    approximate ratios, and the overall units drop out of the g_Gurzhi and g_5
    functions.  We therefore return dimensionless ratios that feed into the
    transport formulae, plus absolute scales where required for plotting.

    Units: densities in 1/nm^3, energies in meV, temperature T in Kelvin.
    """
    mu = p.mu                       # meV
    Tm = p.T * KB_meV_per_K         # temperature in meV
    vF = p.vF                       # nm/ps
    hbar = HBAR_meV_ps              # meV*ps

    # Both nodes (factor 2)
    prefac = 2.0 / (6.0 * np.pi**2 * (hbar * vF)**3)  # 1/(meV^3 * nm^3) per (meV^3)

    # For mu > 0 and finite T, use full result up to T^3 corrections
    n   = prefac * (mu**3 + np.pi**2 * mu * Tm**2)
    # entropy density in 1/(nm^3) with units T/hbar
    s_over_T = (2.0 / (6.0 * np.pi**2 * (hbar * vF)**3)) * (2.0 * np.pi**2 * mu**2 + 14.0/15.0 * np.pi**4 * Tm**2)
    s = s_over_T * Tm  # 1/nm^3 (dimensionless entropy density in natural units)

    # Energy density
    eps = (2.0 / (8.0 * np.pi**2 * (hbar * vF)**3)) * (mu**4 + 2.0*np.pi**2 * mu**2 * Tm**2 + (7.0/15.0) * np.pi**4 * Tm**4)
    P   = eps / 3.0
    w   = eps + P    # w = 4/3 eps (relativistic)

    # Chiral susceptibility (per-node density of states at Fermi level times 2)
    # chi_5 = dn_5/dmu_5 at mu_5=0 = (mu^2 + pi^2 T^2 /3) / (pi^2 (hbar vF)^3)
    chi_5 = (mu**2 + (np.pi**2 / 3.0) * Tm**2) / (np.pi**2 * (hbar * vF)**3)

    return dict(n=n, s=s, eps=eps, P=P, w=w, chi_5=chi_5, Tm=Tm, mu=mu, vF=vF, hbar=hbar)


def hydrodynamic_scales(p: WeylParams):
    """
    Characteristic length scales.  Returns lengths in nm.
    - ell_G  : Gurzhi length  sqrt(eta * vF^2 * tau_mr / w)
    - ell_5  : chiral diffusion length  sqrt(D_5 * tau_5) = vF * sqrt(tau_ee * tau_5 / 3)

    We estimate eta ~ (1/15) * w * tau_ee / vF^2 (standard relativistic Weyl viscosity
    in the hydrodynamic regime; exact prefactor varies by model, e.g. Gorbar et al.).
    For benchmarking purposes we keep this simple.
    """
    thermo = thermodynamic_quantities(p)
    w  = thermo['w']
    vF = p.vF
    tau_mr = p.tau_mr
    tau_ee = p.tau_ee
    tau_5  = p.tau_5

    # Shear viscosity (meV * ps / nm^3)
    eta = (1.0/15.0) * w * tau_ee   # in units of w*tau, length-squared included through vF
    # With w in 1/nm^3 * meV, and tau_ee in ps -> eta in meV ps / nm^3.

    # Gurzhi length
    ell_G = vF * np.sqrt(tau_ee * tau_mr / 15.0)   # = sqrt(eta vF^2 tau / w) with eta = w tau_ee/15 vF^0 form
    # Note: sqrt((w tau_ee/15) * vF^2 * tau_mr / w) = vF * sqrt(tau_ee tau_mr /15)

    # Chiral diffusion length: D_5 = vF^2 tau_ee / 3 in the ballistic-limited case
    ell_5 = vF * np.sqrt(tau_ee * tau_5 / 3.0)

    return dict(eta=eta, ell_G=ell_G, ell_5=ell_5)


def g_gurzhi(x):
    """
    The universal Gurzhi correction factor:
        g(x) = 1 - tanh(x)/x
    where x = W / (2 * ell).

    Limits:
      x -> 0  : g(x) ~ x^2 / 3        (sigma ~ W^2 / (12 ell^2) times bulk)
      x -> inf: g(x) ~ 1 - 1/x        (approaches bulk)
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)

    # Use Taylor for very small x to avoid 0/0
    small = np.abs(x) < 1e-4
    out[small]  = (x[small]**2) / 3.0 - (x[small]**4) * 2.0/15.0
    out[~small] = 1.0 - np.tanh(x[~small]) / x[~small]
    return out


def gurzhi_channel_factor(W_nm, ell_nm):
    """
    g_Gurzhi(W/(2 ell)) : fractional average of Poiseuille flow that carries
    current.  Returns a number in [0, 1].
    """
    return g_gurzhi(W_nm / (2.0 * ell_nm))


# -----------------------------------------------------------------------------
# Quick self-test
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    p = WeylParams()
    th = thermodynamic_quantities(p)
    sc = hydrodynamic_scales(p)

    print(f"mu = {p.mu} meV, T = {p.T} K -> Tm = {th['Tm']:.3f} meV")
    print(f"n = {th['n']:.4e} /nm^3")
    print(f"s = {th['s']:.4e} /nm^3")
    print(f"w = {th['w']:.4e} meV/nm^3")
    print(f"chi_5 = {th['chi_5']:.4e} /(meV nm^3)")
    print()
    print(f"eta = {sc['eta']:.4e} meV ps / nm^3")
    print(f"ell_G = {sc['ell_G']:.3f} nm")
    print(f"ell_5 = {sc['ell_5']:.3f} nm")
    print()
    # Check Gurzhi limits
    for W in [10.0, 50.0, 200.0, 500.0, 2000.0, 10000.0]:
        g = gurzhi_channel_factor(W, sc['ell_G'])
        print(f"W={W:7.1f} nm -> g_Gurzhi(W/2l_G) = {g:.5f}")
