"""
weyl_hydro/transport.py

Hydrodynamic transport coefficients of a Weyl semimetal slab.

We consider a slab of thickness W along y, infinite in x and z.  The driving
fields (uniform in y) are:
   E_x         electric field along x (transport direction)
   dT_x = dT/dx   temperature gradient along x
   B           magnetic field along z (also the chiral-shift axis)
   b_shift     AHE strength (node separation along z)

The relevant equations, in the 'minimal hydrodynamic' regime where the
collision-dominated electron fluid has a single velocity u_x(y), are:

(1) Valley-even Navier-Stokes (momentum balance, steady state, along x):
    eta u_x''(y) - (w / v_F^2 / tau_mr) u_x(y) = F_drive
    where F_drive = -(n e E_x + s dT_x + F_anom)   (driving force density)

    With no-slip BC u_x(+-W/2) = 0:
      u_x(y) = (F_drive / alpha) * [1 - cosh(y/ell_G) / cosh(W/(2 ell_G))]
      where alpha = -(w / v_F^2 / tau_mr) and ell_G = sqrt(eta v_F^2 tau_mr / w).

    Averaged:
      <u_x> = -(F_drive / alpha) * g_Gurzhi(W/(2 ell_G))

(2) Chiral continuity with anomaly (steady state, along x):
    -D_5 mu_5''(y) + mu_5(y)/tau_5 = S_chiral(x-dependences)
    where S_chiral = (e^2 / (2 pi^2 hbar^2 chi_5)) * B * (E_x + beta_grav dT_x)

    With BC mu_5(+-W/2) = 0:
      <mu_5> = S_chiral tau_5 * g_5(W/(2 ell_5))
      ell_5 = sqrt(D_5 tau_5)

(3) Current along x:
    <J_x> = n e <u_x> + sigma_Q (E_x - T dT_x/T)   [simplified intrinsic term]
          + sigma_CME <mu_5> B + sigma_xy^AHE E_y

(4) Heat current along x (open-circuit expression handled later):
    <Q_x> = T s <u_x> + (anomaly- and AHE-induced pieces)

We compute longitudinal thermoelectric tensor elements in the form:

   sigma_xx(W) = d<J_x>/dE_x  at  dT_x = 0
   alpha_xx(W) = -d<J_x>/d(dT_x)  at  E_x = 0            [thermoelectric]
   kappa_xx(W) = -d<Q_x>/d(dT_x)  at  J_x = 0  (open-circuit)
   S_xx(W)     = -(E_x / dT_x)  at  J_x = 0  =  alpha_xx / sigma_xx
   L(W)        = kappa_xx(W) / (T sigma_xx(W))           [Lorenz number]

Benchmarks:
   - B = 0, b = 0, dT_x = 0 :   sigma_xx(W) = sigma_bulk * g_Gurzhi(W/2 ell_G)
                                (standard Gurzhi).
   - W -> infinity with B != 0: recovers bulk hydrodynamic thermoelectrics
                                of Messica--Ostrovsky--Gutman (2023) and
                                Lucas--Davison--Sachdev (2016).
   - dT_x = 0, B != 0, b = 0 :  sigma_xx(W, B) reduces to the anomalous Gurzhi
                                of Sukhachov--Trauzettel (2022): the
                                magnetic-field part of sigma scales with
                                g_5(W/2 ell_5).
"""
from __future__ import annotations
import numpy as np
from params import (
    WeylParams, thermodynamic_quantities, hydrodynamic_scales,
    gurzhi_channel_factor, HBAR_meV_ps, KB_meV_per_K
)


def _coupling_coefficients(p: WeylParams):
    """
    Derived couplings used in the linearized transport equations.

    Returns a dict with:
      n, s, w, chi_5, Tm, vF, hbar              (thermodynamics)
      eta, ell_G, ell_5                          (scales)
      tau_mr, tau_ee, tau_5, W                   (bookkeeping)
      anomaly_coef: eta_A = e^2 / (2 pi^2 hbar^2 chi_5)  (CME-like coefficient)
      beta_grav    : mixed axial-gravitational anomaly coefficient.

    We express 'e' as 1 in natural units since we are reporting transport
    coefficients in arbitrary-but-self-consistent units.  What matters for
    the paper is *relative* behaviour vs W, T, B.
    """
    th  = thermodynamic_quantities(p)
    sc  = hydrodynamic_scales(p)
    out = dict(th)
    out.update(sc)
    out.update(dict(tau_mr=p.tau_mr, tau_ee=p.tau_ee, tau_5=p.tau_5, W=p.W, B=p.B, b=p.b_shift))

    mu = th['mu']
    Tm = th['Tm']

    # Anomaly-induced coefficient relating (E_x B) to the chiral-charge source.
    # The microscopic chiral anomaly gives a source rate
    #    d n_5 / dt|_source = (e^2 / (2 pi^2 hbar^2)) E.B
    # For numerical convenience we work in "natural" units where B is DIMENSIONLESS:
    # B = B_phys / B_*, where B_* = mu^2 / (e hbar vF^2) is the characteristic field
    # at which the cyclotron energy equals mu.  This way B = 1 is "strong".
    # Then the anomaly coefficients become (in our natural units):
    #    kappa_CA  = mu^4 / (2 pi^2 hbar^2 vF^4) -> we absorb constants and track B^2.
    # Equivalently, we express sigma_anom in units where it matches the
    # hydrodynamic scale for reference.
    # For the Sukhachov-Trauzettel-like benchmarks the overall prefactor is an
    # overall scale that cancels in ratios, so we pick a normalisation that
    # puts sigma_anom(B=1) at the same order as sigma_hydro at the reference
    # parameters.
    mu = th['mu']
    vF = p.vF
    # Dimensionless form: numerator has mu^4 / (hbar^2 vF^4) to carry units of charge density squared.
    # This choice ensures sigma_anom^bulk(B=1) ~ n^2 vF^2 tau / w when parameters are typical.
    kappa_CA  = mu**2 / (2.0 * np.pi**2 * HBAR_meV_ps**2 * vF**2)   # source coefficient
    kappa_CME = mu**2 / (4.0 * np.pi**2 * HBAR_meV_ps**2 * vF**2)   # CME current coefficient
    out['kappa_CA']  = kappa_CA
    out['kappa_CME'] = kappa_CME
    # Combined bulk anomalous conductivity: sigma_anom|_bulk = kappa_CME * kappa_CA / chi_5 * B^2 * tau_5
    out['anomaly_coef'] = kappa_CME * kappa_CA / out['chi_5']

    # Mixed axial-gravitational anomaly coefficient.
    # beta_grav = (pi^2 T^2 / 3 mu) contribution at finite mu, plus a
    # mu-independent piece proportional to T^2. In the spirit of
    # Lucas-Davison-Sachdev 2016 and Gooth et al. 2017, the chiral-charge
    # source under a temperature gradient is
    #    ~ (e B / (2 pi^2 hbar^2 chi_5)) * beta_grav * dT_x
    # We take the combined form:
    #    beta_grav = (pi^2 Tm^2) / (3 mu + epsilon)   for finite mu
    # At the neutrality point mu -> 0, we use a finite-T form beta_grav ~ (pi/3) Tm.
    if mu > 1e-6 * Tm:
        beta_grav = (np.pi**2 * Tm**2) / (3.0 * mu + 1e-9)
    else:
        beta_grav = (np.pi / 3.0) * Tm
    out['beta_grav'] = beta_grav

    # Anomalous Hall coefficient: sigma_xy^AHE in natural units.
    # sigma_xy^AHE = e^2 b / (2 pi^2 hbar) * (arbitrary overall unit).
    # We set e = 1 and carry it as an overall normalization.
    out['sigma_ahe'] = p.b_shift / (2.0 * np.pi**2 * HBAR_meV_ps)

    # Intrinsic ('quantum-critical' / incoherent) conductivity sigma_Q.
    # For a clean Weyl fluid, sigma_Q ~ (e^2 / hbar) * (T/(hbar vF)) with a
    # numerical coefficient.  We use a standard O(1) prefactor.
    out['sigma_Q'] = (Tm / (HBAR_meV_ps * p.vF)) / (6.0 * np.pi)

    return out


def _drude_like_prefactors(coeffs):
    """
    Return the 'force -> velocity' conversion factors and their bulk limits.

    The Navier-Stokes equation
       eta u'' - (w / v_F^2 tau_mr) u = F
    has bulk solution (W -> infinity): u = -F v_F^2 tau_mr / w.
    In a slab with no-slip BC:
       <u_x> = -F v_F^2 tau_mr / w * g_Gurzhi(W/2 ell_G).
    """
    w  = coeffs['w']
    vF = coeffs['vF']
    tau_mr = coeffs['tau_mr']
    ell_G  = coeffs['ell_G']
    W      = coeffs['W']

    g = gurzhi_channel_factor(W, ell_G)  # 0 <= g <= 1
    # <u_x> = - F / alpha * g,  with alpha = w/(v_F^2 tau_mr)
    u_over_F_bulk = (vF * vF * tau_mr) / w   # velocity response per unit (outward) force, bulk
    # Bulk: g -> 1 when W >> ell_G.  For finite W: multiply by g.
    return u_over_F_bulk, g


def _chiral_suppression(coeffs):
    """
    Return the chiral-charge suppression factor g_5(W/2 ell_5) and the bulk
    response.

    Chiral continuity equation (steady state, along x with uniform drivers):
        mu_5/tau_5 - D_5 mu_5''(y) = S_chiral
    with BC mu_5(+-W/2) = 0:
        <mu_5> = S_chiral tau_5 * g_5(W/(2 ell_5))
    """
    ell_5 = coeffs['ell_5']
    W     = coeffs['W']
    g5 = gurzhi_channel_factor(W, ell_5)
    return coeffs['tau_5'], g5


def bulk_conductivity(p: WeylParams):
    """
    Bulk (W -> inf) transport coefficients, for benchmarking against
    Messica et al. / Lucas-Davison-Sachdev / Boltzmann Drude.
    """
    c = _coupling_coefficients(p)
    n = c['n']; w = c['w']; vF = c['vF']; tau_mr = c['tau_mr']

    # Drude-like hydro: sigma_Drude = n^2 e^2 v_F^2 tau_mr / w
    sigma_drude = (n**2) * (vF**2) * tau_mr / w

    # Anomalous (CME-induced) bulk conductivity: proportional to B^2
    sigma_anom_per_B2 = c['kappa_CME'] * c['kappa_CA'] * c['tau_5'] / c['chi_5']

    # Thermopower (bulk): S = -(alpha / sigma). For the hydrodynamic Weyl gas,
    # alpha_bulk ~ s * vF^2 tau_mr / w * n  (from Navier-Stokes with dT drive),
    # which simplifies after dividing by sigma_drude:
    #   S_bulk = -(alpha/sigma) = -s / (n e)     (Wiedemann Franz-respecting form at low T)
    # In practice finite corrections arise from sigma_Q and the anomaly.
    # We keep the leading piece here:
    S_bulk = -c['s'] / (n + 1e-30)

    # Thermal conductivity bulk:
    #   kappa_bulk (closed circuit) = s^2 T vF^2 tau_mr / w  +  other pieces
    kappa_bulk = (c['s']**2) * c['Tm'] * (vF**2) * tau_mr / w

    # Lorenz number
    L_bulk = kappa_bulk / (c['Tm'] * sigma_drude + 1e-30)

    return dict(
        sigma_bulk=sigma_drude,
        sigma_anom_per_B2=sigma_anom_per_B2,
        S_bulk=S_bulk,
        kappa_bulk=kappa_bulk,
        L_bulk=L_bulk,
        coeffs=c,
    )


def channel_transport(p: WeylParams):
    """
    Full slab-geometry transport coefficients.  This is the main function.

    We compute:
       sigma_xx(W)   : longitudinal electrical conductivity.
       alpha_xx(W)   : longitudinal thermoelectric conductivity (  d<J_x>/d(-dT_x)  )
       kappa_xx(W)   : longitudinal thermal conductivity (open circuit).
       sigma_ahe_eff : effective AHE (geometry-modulated).
       alpha_xy(W)   : transverse Nernst coefficient.
       S_xx(W)       : Seebeck coefficient  = alpha_xx / sigma_xx.
       L(W)          : Lorenz ratio  = kappa_xx / (T sigma_xx).

    Breakdown:
       sigma_xx(W)   = sigma_hydro(W) + sigma_anom(W, B)
         sigma_hydro(W) = sigma_bulk^hydro  * g_Gurzhi(W/2 ell_G)
         sigma_anom(W)  = sigma_anom_bulk(B) * g_5(W/2 ell_5)
                         ( Sukhachov-Trauzettel anomalous-Gurzhi form )
       alpha_xx(W)   = alpha_hydro(W) + alpha_anom(W, B)
         alpha_hydro(W) = alpha_bulk^hydro * g_Gurzhi(W/2 ell_G)
         alpha_anom(W)  = (anomaly_coef * B)(anomaly_coef * B * beta_grav * tau_5)
                          * g_5(W/2 ell_5)
            == magnitude tracks the mixed axial-gravitational anomaly piece.
       kappa_xx(W)   has an 'open-circuit' piece that cancels with sigma_xx.

    Physics subtlety:  sigma_hydro and sigma_anom have *different* spatial
    profiles (viscous Poiseuille vs chiral-charge diffusion), so they carry
    independent geometric factors g_Gurzhi and g_5.  This is what produces
    the decoupling between longitudinal Gurzhi and the 'chiral Gurzhi'.
    """
    c = _coupling_coefficients(p)
    n = c['n']; w = c['w']; vF = c['vF']; tau_mr = c['tau_mr']
    Tm = c['Tm']; s = c['s']
    B  = p.B

    # Gurzhi factors
    g_G = gurzhi_channel_factor(c['W'], c['ell_G'])
    g_5 = gurzhi_channel_factor(c['W'], c['ell_5'])

    # --- Electrical conductivity ---
    # Hydrodynamic piece (drive = n e E_x on fluid):
    sigma_hydro_bulk = (n**2) * (vF**2) * tau_mr / w
    sigma_hydro = sigma_hydro_bulk * g_G

    # Anomaly piece (chiral-anomaly-induced CME, scales as B^2).
    # sigma_anom = kappa_CME * <mu_5> * B / E_x
    #            = kappa_CME * (kappa_CA B E_x tau_5 / chi_5) * B / E_x * g_5
    #            = kappa_CME * kappa_CA * B^2 * tau_5 / chi_5 * g_5
    sigma_anom_bulk = c['kappa_CME'] * c['kappa_CA'] * (B**2) * c['tau_5'] / c['chi_5']
    sigma_anom = sigma_anom_bulk * g_5

    # Plus an incoherent 'sigma_Q' background that does *not* couple to u_x
    sigma_xx = sigma_hydro + sigma_anom + c['sigma_Q']

    # --- Thermoelectric (alpha_xx) ---
    # Hydrodynamic piece (drive = s (-dT_x) on fluid):
    # from Navier-Stokes with drive F = -s dT_x, current = n e <u_x>:
    alpha_hydro_bulk = n * s * (vF**2) * tau_mr / w
    alpha_hydro = alpha_hydro_bulk * g_G

    # Mixed axial-gravitational anomaly piece:
    # chiral source from (B * dT_x) gives extra CME current,
    # with the beta_grav factor encoding the axial-gravitational contribution.
    # alpha_anom_bulk = kappa_CME * kappa_CA * B^2 * tau_5 * beta_grav / chi_5
    alpha_anom_bulk = c['kappa_CME'] * c['kappa_CA'] * (B**2) * c['tau_5'] * c['beta_grav'] / c['chi_5']
    alpha_anom = alpha_anom_bulk * g_5

    alpha_xx = alpha_hydro + alpha_anom

    # --- Thermal conductivity (open circuit) ---
    # In the hydrodynamic regime: Q_x = T s u_x (leading).  So kappa_closed = T s^2 vF^2 tau_mr / w * g_G.
    # Plus an anomaly-induced kappa_anom proportional to B^2 * T * beta_grav^2 * tau_5 * g_5.
    kappa_closed_hydro = Tm * (s**2) * (vF**2) * tau_mr / w * g_G
    kappa_closed_anom  = (Tm * c['kappa_CME'] * c['kappa_CA']
                          * (B**2) * (c['beta_grav']**2) * c['tau_5'] / c['chi_5'] * g_5)
    kappa_closed = kappa_closed_hydro + kappa_closed_anom

    # Open circuit (thermal): kappa_oc = kappa_closed - Tm * alpha_xx^2 / sigma_xx
    kappa_xx = kappa_closed - Tm * (alpha_xx**2) / (sigma_xx + 1e-30)

    # --- Seebeck and Lorenz ---
    S_xx = alpha_xx / (sigma_xx + 1e-30)
    L    = kappa_xx / (Tm * sigma_xx + 1e-30)

    # --- Transverse Nernst from AHE times hydrodynamic flow pattern ---
    # The key non-graphene piece: sigma_xy^AHE multiplies the 'flow-induced' voltage.
    # In the slab with uniform fields, alpha_xy picks up a (-sigma_ahe * S_xx) kind of term,
    # whose W-dependence tracks g_G through S_xx.
    alpha_xy = c['sigma_ahe'] * S_xx    # simple schematic, sensitive to g_G via S_xx

    return dict(
        sigma_xx=sigma_xx,
        sigma_hydro=sigma_hydro,
        sigma_anom=sigma_anom,
        sigma_Q=c['sigma_Q'],
        alpha_xx=alpha_xx,
        alpha_hydro=alpha_hydro,
        alpha_anom=alpha_anom,
        kappa_xx=kappa_xx,
        kappa_closed=kappa_closed,
        S_xx=S_xx,
        L=L,
        alpha_xy=alpha_xy,
        g_G=g_G,
        g_5=g_5,
        ell_G=c['ell_G'],
        ell_5=c['ell_5'],
    )


def seebeck_bulk_reference(p: WeylParams):
    """
    The bulk Seebeck reference value (W -> infinity), used for normalisation
    in plots.  Obtained by replacing g_G = g_5 = 1.
    """
    p2 = WeylParams(**{**p.__dict__, 'W': 1.0e12})  # 1e12 nm effectively bulk
    out = channel_transport(p2)
    return out['S_xx'], out['sigma_xx'], out['alpha_xx'], out['kappa_xx'], out['L']


# -----------------------------------------------------------------------------
# Quick self-consistency checks
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Benchmark 1: Gurzhi at B = 0, b = 0.
    print("=== Benchmark 1: Gurzhi (B=0, b=0) ===")
    p0 = WeylParams(B=0.0, b_shift=0.0)
    bulk = bulk_conductivity(p0)
    print(f"  sigma_bulk = {bulk['sigma_bulk']:.4e}")
    print(f"  S_bulk     = {bulk['S_bulk']:.4e}")

    for W in [10, 50, 200, 500, 2000, 1e4]:
        p0.W = W
        out = channel_transport(p0)
        print(f"  W={W:7.1f} nm:  sigma/sigma_bulk = {out['sigma_hydro']/bulk['sigma_bulk']:.4f} "
              f"(should equal g_G = {out['g_G']:.4f})")

    # Benchmark 2: Bulk Seebeck recovers S = s/n form.
    print("\n=== Benchmark 2: Bulk W -> inf Seebeck ===")
    p1 = WeylParams(B=0.0, b_shift=0.0, W=1.0e10)
    out = channel_transport(p1)
    th = thermodynamic_quantities(p1)
    print(f"  S_bulk from channel_transport = {out['S_xx']:.4e}")
    print(f"  analytic  s/n                 = {th['s']/th['n']:.4e}")
    print(f"  ratio     = {out['S_xx'] / (th['s']/th['n']):.4f} "
          f"(should be ~1, small deviation from sigma_Q background)")

    # Benchmark 3: Anomalous Gurzhi - magnetic part g_5 dependence on W.
    print("\n=== Benchmark 3: Anomalous Gurzhi sigma_anom(W)/sigma_anom_bulk = g_5 ===")
    p2 = WeylParams(B=1.0, b_shift=0.0)
    for W in [10, 50, 200, 1000, 5000, 1e4]:
        p2.W = W
        out = channel_transport(p2)
        ratio = out['sigma_anom'] / ((_coupling_coefficients(p2)['anomaly_coef']**2) * (p2.B**2) * p2.tau_5 + 1e-30)
        print(f"  W={W:7.1f} nm:  sigma_anom(W)/bulk = {ratio:.4f} (g_5 = {out['g_5']:.4f})")
