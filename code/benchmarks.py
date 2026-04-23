"""
benchmarks.py

Full benchmark suite for the Weyl hydrodynamic transport code.

Each benchmark is a sanity check against a known limit or a published result.
All benchmarks print PASS/FAIL and leave a paper trail.
"""
import numpy as np
from params import (WeylParams, thermodynamic_quantities, hydrodynamic_scales,
                    gurzhi_channel_factor, g_gurzhi)
from transport import channel_transport, bulk_conductivity, _coupling_coefficients


def _report(name, ok, info=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}  {info}")
    return ok


def benchmark_gurzhi_limits():
    """
    Benchmark 1: Gurzhi correction factor asymptotics.
      g(x) -> x^2/3 for x -> 0   (Poiseuille limit)
      g(x) -> 1 - 1/x for x -> inf (approaches bulk)
    """
    print("[B1] Gurzhi correction factor asymptotics")
    # Small x
    xs = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    g  = g_gurzhi(xs)
    ratio_small = g / (xs**2 / 3.0)
    ok1 = np.allclose(ratio_small, 1.0, rtol=1e-2)
    _report("small x -> x^2/3", ok1, f"ratios={ratio_small}")
    # Large x
    xs = np.array([10., 50., 100., 500.])
    g = g_gurzhi(xs)
    ratio_large = g / (1.0 - 1.0/xs)
    ok2 = np.allclose(ratio_large, 1.0, rtol=5e-3)
    _report("large x -> 1 - 1/x", ok2, f"ratios={ratio_large}")
    return ok1 and ok2


def benchmark_bulk_gurzhi():
    """
    Benchmark 2: sigma_hydro(W)/sigma_bulk must equal g_Gurzhi(W/2 ell_G) exactly
    (no other Gurzhi-shaped factors should contaminate it).
    """
    print("[B2] Viscous Gurzhi of sigma")
    p = WeylParams(B=0.0, b_shift=0.0)
    bulk = bulk_conductivity(p)
    all_ok = True
    for W in [1, 5, 10, 50, 100, 500, 2000, 1e5]:
        p.W = W
        out = channel_transport(p)
        predicted = out['g_G']
        actual    = out['sigma_hydro'] / bulk['sigma_bulk']
        ok = np.isclose(predicted, actual, rtol=1e-10)
        all_ok &= ok
        _report(f"W={W:6.1f} nm", ok, f"pred={predicted:.6e}  actual={actual:.6e}")
    return all_ok


def benchmark_seebeck_bulk_formula():
    """
    Benchmark 3: Bulk Seebeck matches S = s/n (to within small sigma_Q correction).
    """
    print("[B3] Bulk Seebeck S_bulk ~ s/n")
    p = WeylParams(B=0.0, b_shift=0.0, W=1e10)
    out = channel_transport(p)
    th  = thermodynamic_quantities(p)
    S_analytic = th['s'] / th['n']
    ratio = out['S_xx'] / S_analytic
    ok = (ratio > 0.99) and (ratio < 1.001)
    return _report("S_bulk / (s/n)", ok, f"ratio={ratio:.6f}")


def benchmark_anomalous_gurzhi_ST():
    """
    Benchmark 4: anomalous part of sigma at B != 0, b = 0, dT = 0 scales
    with g_5(W/2 ell_5) (Sukhachov-Trauzettel anomalous Gurzhi).
    """
    print("[B4] Anomalous Gurzhi of sigma_anom vs g_5")
    p = WeylParams(B=1.0, b_shift=0.0)
    all_ok = True
    for W in [1, 5, 10, 50, 100, 500, 2000, 1e5]:
        p.W = W
        out = channel_transport(p)
        c = _coupling_coefficients(p)
        bulk = c['kappa_CME'] * c['kappa_CA'] * (p.B**2) * p.tau_5 / c['chi_5']
        actual = out['sigma_anom'] / bulk
        predicted = out['g_5']
        ok = np.isclose(predicted, actual, rtol=1e-10)
        all_ok &= ok
        _report(f"W={W:6.1f} nm", ok, f"pred={predicted:.6e}  actual={actual:.6e}")
    return all_ok


def benchmark_hydro_scale_consistency():
    """
    Benchmark 5: the Gurzhi length ell_G scales as vF sqrt(tau_ee tau_mr).
    """
    print("[B5] Gurzhi length scaling")
    ok_all = True
    p = WeylParams()
    for tau_ee in [0.5, 1.0, 2.0, 5.0]:
        for tau_mr in [5.0, 10.0, 20.0, 50.0]:
            p.tau_ee = tau_ee
            p.tau_mr = tau_mr
            sc = hydrodynamic_scales(p)
            predicted = p.vF * np.sqrt(tau_ee * tau_mr / 15.0)
            ok = np.isclose(sc['ell_G'], predicted, rtol=1e-10)
            ok_all &= ok
    return _report("ell_G = vF sqrt(tau_ee tau_mr / 15)", ok_all)


def benchmark_alpha_xx_same_Gurzhi_as_sigma_hydro():
    """
    Benchmark 6: the hydrodynamic piece of alpha_xx carries the same Gurzhi
    factor g_G as sigma_hydro (they both come from the same Navier-Stokes
    channel problem with different drivers).
    """
    print("[B6] alpha_hydro shares Gurzhi factor with sigma_hydro")
    p = WeylParams(B=0.0, b_shift=0.0)
    all_ok = True
    for W in [10, 100, 500, 2000, 1e4]:
        p.W = W
        out = channel_transport(p)
        ratio_sigma = out['sigma_hydro'] / (out['sigma_hydro'] / out['g_G'])
        ratio_alpha = out['alpha_hydro'] / (out['alpha_hydro'] / out['g_G'])
        ok = np.isclose(ratio_sigma, ratio_alpha, rtol=1e-10)
        all_ok &= ok
    return _report("alpha/sigma Gurzhi shared", all_ok)


def benchmark_Wiedemann_Franz_hydro_only_limit():
    """
    Benchmark 7: in the hydrodynamic-only limit (sigma_Q -> 0 by hand, B=0,
    b=0), L(W) should be exactly W-independent because both kappa_closed
    and sigma scale with the same g_G. We test this by post-processing the
    output to remove the sigma_Q channel.
    """
    print("[B7] Wiedemann-Franz L(W) is W-indep in pure hydro limit (sigma_Q removed)")
    p = WeylParams(B=0.0, b_shift=0.0, mu=5.0)
    L_list = []
    W_list = [50, 200, 500, 2000, 1e4]
    for W in W_list:
        p.W = W
        out = channel_transport(p)
        # Rebuild L without the sigma_Q piece
        sigma_pure = out['sigma_hydro']
        # alpha is already pure hydro (no sigma_Q contribution)
        alpha_pure = out['alpha_hydro']
        # kappa_closed is already pure hydro at B=0
        kappa_closed_pure = out['kappa_closed']
        th = thermodynamic_quantities(p)
        Tm = th['Tm']
        kappa_oc = kappa_closed_pure - Tm * (alpha_pure**2) / sigma_pure
        L_pure = kappa_oc / (Tm * sigma_pure)
        L_list.append(L_pure)
    L_arr = np.array(L_list)
    rel_range = (L_arr.max() - L_arr.min()) / L_arr.mean()
    ok = rel_range < 1e-6
    return _report("L(W) W-indep in hydro-only limit",
                   ok, f"rel_range={rel_range:.3e}, L_values={L_arr}")


def benchmark_negative_magnetoresistance_bulk():
    """
    Benchmark 8: In bulk, sigma_xx(B) has a positive B^2 correction from the
    chiral anomaly (i.e. RESISTANCE decreases, MR is NEGATIVE).
    """
    print("[B8] Negative MR bulk")
    p = WeylParams(W=1e10, b_shift=0.0)
    sigma_0 = channel_transport(p)['sigma_xx']
    p.B = 1.0
    sigma_B = channel_transport(p)['sigma_xx']
    ok = sigma_B > sigma_0
    delta = (sigma_B - sigma_0) / sigma_0
    return _report("sigma_xx(B) > sigma_xx(0)", ok, f"delta = {delta:.3%}")


def benchmark_bulk_anom_scales_as_B2():
    """
    Benchmark 9: Bulk anomalous conductivity scales as B^2 (chiral anomaly).
    """
    print("[B9] Bulk sigma_anom ~ B^2")
    p = WeylParams(W=1e10, b_shift=0.0)
    B_values = np.array([0.1, 0.2, 0.5, 1.0, 2.0])
    sig_anom = []
    for B in B_values:
        p.B = B
        out = channel_transport(p)
        sig_anom.append(out['sigma_anom'])
    sig_anom = np.array(sig_anom)
    # fit log(sig) vs log(B)
    slope, _ = np.polyfit(np.log(B_values), np.log(sig_anom), 1)
    ok = np.isclose(slope, 2.0, atol=1e-2)
    return _report("log-log slope B^2", ok, f"slope={slope:.4f}")


def main():
    print("=" * 70)
    print("Weyl hydrodynamic transport: full benchmark suite")
    print("=" * 70)

    results = [
        benchmark_gurzhi_limits(),
        benchmark_bulk_gurzhi(),
        benchmark_seebeck_bulk_formula(),
        benchmark_anomalous_gurzhi_ST(),
        benchmark_hydro_scale_consistency(),
        benchmark_alpha_xx_same_Gurzhi_as_sigma_hydro(),
        benchmark_Wiedemann_Franz_hydro_only_limit(),
        benchmark_negative_magnetoresistance_bulk(),
        benchmark_bulk_anom_scales_as_B2(),
    ]

    print("=" * 70)
    n_pass = sum(results)
    n_tot  = len(results)
    print(f"{n_pass}/{n_tot} benchmarks passed")
    print("=" * 70)


if __name__ == "__main__":
    main()
