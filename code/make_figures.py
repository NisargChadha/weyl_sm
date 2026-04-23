"""
make_figures.py

Generate all figures for the paper "Thermoelectric Gurzhi effect in Weyl semimetals".
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from params import WeylParams, thermodynamic_quantities, hydrodynamic_scales, gurzhi_channel_factor, g_gurzhi
from transport import channel_transport, _coupling_coefficients, bulk_conductivity

# Style
plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 250,
    "font.size": 10,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": False,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "figure.constrained_layout.use": True,
})

FIG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figures"))
os.makedirs(FIG_DIR, exist_ok=True)


def save(name):
    path = os.path.join(FIG_DIR, name)
    plt.savefig(path, bbox_inches="tight")
    print(f"  wrote {path}")


# =============================================================================
# Figure 1: the Gurzhi correction factor g(x), benchmarks and schematic
# =============================================================================
def figure_1_gurzhi_factor():
    """Figure 1: the universal Gurzhi function g(x) = 1 - tanh(x)/x."""
    fig, ax = plt.subplots(1, 1, figsize=(4.2, 3.2))
    x = np.logspace(-2, 2, 400)
    g = g_gurzhi(x)
    ax.loglog(x, g, 'k-', lw=2.0, label=r"$g(x) = 1 - \tanh(x)/x$")
    ax.loglog(x, x**2/3.0, 'r--', lw=1.0, label=r"$x^2/3$ (Poiseuille)")
    ax.loglog(x, np.ones_like(x), 'b--', lw=1.0, alpha=0.6, label="bulk limit")
    ax.set_xlabel(r"$x = W / 2\ell$")
    ax.set_ylabel(r"$g(x)$")
    ax.set_title("Gurzhi correction factor")
    ax.legend(loc='lower right')
    ax.grid(True, which='both', alpha=0.3)
    save("fig1_gurzhi_factor.png")
    plt.close()


# =============================================================================
# Figure 2: benchmark - sigma_hydro / sigma_bulk = g_Gurzhi(W/2 ell_G)
# and sigma_anom(W) / sigma_anom^bulk = g_5(W/2 ell_5)
# =============================================================================
def figure_2_viscous_and_chiral_gurzhi():
    """Figure 2: two independent Gurzhi scales in a WSM slab."""
    p = WeylParams(B=0.2, b_shift=0.0)
    W_values = np.logspace(0, 4.5, 200)   # 1 to 30,000 nm

    sig_hydro, sig_anom, g_G, g_5 = [], [], [], []
    for W in W_values:
        p.W = W
        out = channel_transport(p)
        sig_hydro.append(out['sigma_hydro'])
        sig_anom.append(out['sigma_anom'])
        g_G.append(out['g_G'])
        g_5.append(out['g_5'])

    sig_hydro = np.array(sig_hydro)
    sig_anom  = np.array(sig_anom)
    g_G = np.array(g_G); g_5 = np.array(g_5)

    # Normalize
    sig_hydro_bulk = sig_hydro[-1]
    sig_anom_bulk  = sig_anom[-1]

    fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.5))
    ax.semilogx(W_values, sig_hydro / sig_hydro_bulk, 'b-', lw=2.0,
                label=r"$\sigma_{\rm hydro}(W) / \sigma_{\rm hydro}^\infty$")
    ax.semilogx(W_values, sig_anom  / sig_anom_bulk,  'r-', lw=2.0,
                label=r"$\sigma_{\rm anom}(W) / \sigma_{\rm anom}^\infty$")
    # Analytical reference curves
    c = _coupling_coefficients(p)
    ax.semilogx(W_values, g_gurzhi(W_values / (2 * c['ell_G'])), 'b--', lw=1.0, alpha=0.6,
                label=r"$g(W/2\ell_G)$")
    ax.semilogx(W_values, g_gurzhi(W_values / (2 * c['ell_5'])), 'r--', lw=1.0, alpha=0.6,
                label=r"$g(W/2\ell_5)$")

    ax.axvline(2 * c['ell_G'], color='b', alpha=0.3, lw=0.8)
    ax.axvline(2 * c['ell_5'], color='r', alpha=0.3, lw=0.8)
    ax.text(2 * c['ell_G'] * 1.1, 0.05, r"$2\ell_G$", color='b', fontsize=9)
    ax.text(2 * c['ell_5'] * 1.1, 0.05, r"$2\ell_5$", color='r', fontsize=9)

    ax.set_xlabel(r"Slab width $W$ (nm)")
    ax.set_ylabel("Normalized conductivity")
    ax.set_title("Two independent Gurzhi scales")
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(-0.05, 1.05)
    save("fig2_two_gurzhi_scales.png")
    plt.close()


# =============================================================================
# Figure 3: Seebeck and Lorenz number vs W, at B=0, varying mu
# (Thermoelectric Gurzhi: main new result)
# =============================================================================
def figure_3_thermoelectric_gurzhi():
    """Main figure: Seebeck S(W) / S_bulk and Lorenz L(W)/L_bulk vs W."""
    W_values = np.logspace(0.5, 4.2, 200)

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))

    colors = plt.cm.viridis(np.linspace(0.1, 0.8, 4))
    mu_values = [2.0, 5.0, 10.0, 20.0]
    B = 0.0
    for mu, col in zip(mu_values, colors):
        p = WeylParams(B=B, b_shift=0.0, mu=mu)
        S_arr, L_arr = [], []
        for W in W_values:
            p.W = W
            out = channel_transport(p)
            S_arr.append(out['S_xx'])
            L_arr.append(out['L'])
        S_arr = np.array(S_arr); L_arr = np.array(L_arr)
        S_bulk = S_arr[-1]
        axes[0].semilogx(W_values, S_arr/S_bulk, color=col, lw=1.8,
                         label=fr"$\mu={mu:.1f}$ meV")
        axes[1].semilogx(W_values, L_arr, color=col, lw=1.8,
                         label=fr"$\mu={mu:.1f}$ meV")

    axes[0].set_xlabel(r"Slab width $W$ (nm)")
    axes[0].set_ylabel(r"$S_{xx}(W) / S_{xx}^\infty$")
    axes[0].set_title("Thermoelectric Gurzhi: Seebeck suppression")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, which='both', alpha=0.3)

    axes[1].set_xlabel(r"Slab width $W$ (nm)")
    axes[1].set_ylabel(r"$L(W) = \kappa_{xx}/(T\,\sigma_{xx})$")
    axes[1].set_title("Lorenz number in narrow slab")
    axes[1].set_yscale('log')
    axes[1].legend(fontsize=9)
    axes[1].grid(True, which='both', alpha=0.3)
    save("fig3_thermoelectric_gurzhi.png")
    plt.close()


# =============================================================================
# Figure 4: L(W, T) contour plot  -- the "thermoelectric Gurzhi map"
# =============================================================================
def figure_4_L_contour():
    """2D contour of Lorenz number in the (W, T) plane."""
    W_values = np.logspace(1, 4, 70)
    T_values = np.linspace(2, 40, 50)

    mu = 5.0
    L_grid = np.zeros((len(T_values), len(W_values)))
    SG_grid = np.zeros_like(L_grid)
    S_grid = np.zeros_like(L_grid)
    for i, T in enumerate(T_values):
        for j, W in enumerate(W_values):
            p = WeylParams(B=0.0, b_shift=0.0, mu=mu, T=T, W=W)
            out = channel_transport(p)
            L_grid[i, j] = out['L']
            SG_grid[i, j] = out['sigma_hydro'] / (out['sigma_xx'])
            S_grid[i, j] = np.abs(out['S_xx'])

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    X, Y = np.meshgrid(W_values, T_values)

    # L / L0 where L0 = pi^2/3 * (kB/e)^2 (in our units, kB = 1, so L0 = pi^2/3)
    L0 = np.pi**2 / 3.0
    pcm = axes[0].pcolormesh(X, Y, L_grid/L0, cmap='viridis',
                             vmin=0, vmax=np.nanpercentile(L_grid/L0, 98), shading='auto')
    axes[0].set_xscale('log')
    axes[0].set_xlabel(r"$W$ (nm)")
    axes[0].set_ylabel(r"$T$ (K)")
    axes[0].set_title(r"Lorenz ratio $L / L_0$")
    plt.colorbar(pcm, ax=axes[0], label=r"$L/L_0$")

    # Seebeck |S_xx|
    pcm2 = axes[1].pcolormesh(X, Y, S_grid, cmap='magma', shading='auto')
    axes[1].set_xscale('log')
    axes[1].set_xlabel(r"$W$ (nm)")
    axes[1].set_ylabel(r"$T$ (K)")
    axes[1].set_title(r"Seebeck $|S_{xx}|$ (arb. units)")
    plt.colorbar(pcm2, ax=axes[1], label=r"$|S_{xx}|$")

    save("fig4_L_map.png")
    plt.close()


# =============================================================================
# Figure 5: the chiral-anomaly contribution to Seebeck
# Shows alpha_anom/alpha_hydro as a function of W and B
# =============================================================================
def figure_5_anom_contribution():
    """How much does the axial-gravitational anomaly contribute vs hydrodynamic?"""
    W_values = np.logspace(1, 4.5, 150)
    B_values = [0.1, 0.3, 0.5, 1.0]

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(B_values)))

    # --- Left panel: width dependence of the anomaly-induced chiral-anomaly
    # contribution to alpha, normalized to bulk ---
    for B, col in zip(B_values, colors):
        p = WeylParams(B=B, b_shift=0.0, mu=5.0)
        alpha_anom_W = []
        for W in W_values:
            p.W = W
            out = channel_transport(p)
            alpha_anom_W.append(out['alpha_anom'])
        alpha_anom_W = np.array(alpha_anom_W)
        alpha_anom_bulk = alpha_anom_W[-1]
        axes[0].semilogx(W_values, alpha_anom_W / alpha_anom_bulk,
                         color=col, lw=1.8, label=fr"$B={B}$")

    # Reference curve: g_5(W/2 ell_5)
    p0 = WeylParams(B=0.1, b_shift=0.0, mu=5.0)
    sc = hydrodynamic_scales(p0)
    axes[0].semilogx(W_values, g_gurzhi(W_values / (2 * sc['ell_5'])), 'k--', lw=1.2,
                     label=r"$g(W/2\ell_5)$")
    axes[0].axvline(2 * sc['ell_5'], color='gray', lw=0.8, alpha=0.5)
    axes[0].text(2 * sc['ell_5'] * 1.15, 0.05, r"$2\ell_5$", color='gray', fontsize=9)
    axes[0].set_xlabel(r"Slab width $W$ (nm)")
    axes[0].set_ylabel(r"$\alpha_{\rm anom}(W) / \alpha_{\rm anom}^\infty$")
    axes[0].set_title("Mixed axial-gravitational signal\n vs channel width")
    axes[0].legend(fontsize=9, loc='upper left')
    axes[0].grid(True, which='both', alpha=0.3)

    # Right panel: Seebeck vs B at fixed W
    W_fixed = [200.0, 500.0, 1500.0, 5000.0]
    B_range = np.linspace(0, 1.5, 60)
    colors = plt.cm.coolwarm(np.linspace(0.2, 0.85, len(W_fixed)))
    for W, col in zip(W_fixed, colors):
        S_arr = []
        for B in B_range:
            p = WeylParams(B=B, b_shift=0.0, mu=5.0, W=W)
            out = channel_transport(p)
            S_arr.append(out['S_xx'])
        S_arr = np.array(S_arr)
        axes[1].plot(B_range, S_arr / S_arr[0], color=col, lw=1.8, label=fr"$W={W:.0f}$ nm")

    axes[1].set_xlabel(r"$B$ (dim-less)")
    axes[1].set_ylabel(r"$S_{xx}(B) / S_{xx}(0)$")
    axes[1].set_title(r"Anomaly-induced $S_{xx}(B)$" + "\n width-dependent")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, which='both', alpha=0.3)

    save("fig5_anomaly_contribution.png")
    plt.close()


# =============================================================================
# Figure 6: Poiseuille profile u_x(y) analytical
# =============================================================================
def figure_6_poiseuille_profiles():
    """Profile of u_x(y) in the slab for various W/ell_G ratios."""
    fig, ax = plt.subplots(1, 1, figsize=(4.6, 3.3))
    p = WeylParams()
    sc = hydrodynamic_scales(p)
    ell_G = sc['ell_G']

    y_rel = np.linspace(-0.5, 0.5, 200)  # y/W
    Ws = [0.5 * ell_G, 1.0 * ell_G, 2.0 * ell_G, 5.0 * ell_G, 20.0 * ell_G]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(Ws)))
    for W, col in zip(Ws, colors):
        y = y_rel * W
        u = 1.0 - np.cosh(y / ell_G) / np.cosh(W / (2 * ell_G))
        ax.plot(y_rel, u, color=col, lw=1.8, label=fr"$W/\ell_G = {W/ell_G:.1f}$")
    ax.set_xlabel(r"$y/W$")
    ax.set_ylabel(r"$u_x(y)\,/\,u_{\rm max}^{\rm bulk}$")
    ax.set_title("Hydrodynamic flow profile")
    ax.legend(fontsize=9, loc='lower center')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(-0.5, 0.5)
    save("fig6_poiseuille.png")
    plt.close()


# =============================================================================
# Figure 7: Temperature dependence at fixed W (thermoelectric "Gurzhi minimum")
# =============================================================================
def figure_7_T_dependence():
    """Non-monotonic temperature dependence of L and S at fixed W."""
    T_values = np.linspace(2, 50, 80)
    W_values = [200, 500, 1500, 5000]

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))
    colors = plt.cm.coolwarm(np.linspace(0.15, 0.85, len(W_values)))
    for W, col in zip(W_values, colors):
        L_arr, S_arr = [], []
        for T in T_values:
            p = WeylParams(B=0.3, b_shift=0.0, mu=5.0, T=T, W=W)
            out = channel_transport(p)
            L_arr.append(out['L'])
            S_arr.append(out['S_xx'])
        L_arr = np.array(L_arr); S_arr = np.array(S_arr)
        axes[0].plot(T_values, L_arr / (np.pi**2 / 3.0), color=col, lw=1.8, label=fr"$W={W}$ nm")
        axes[1].plot(T_values, S_arr, color=col, lw=1.8, label=fr"$W={W}$ nm")

    axes[0].set_xlabel(r"$T$ (K)")
    axes[0].set_ylabel(r"$L / L_0$")
    axes[0].set_title(r"Lorenz ratio at $B=0.3$")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, which='both', alpha=0.3)
    axes[0].set_yscale('log')

    axes[1].set_xlabel(r"$T$ (K)")
    axes[1].set_ylabel(r"$S_{xx}$ (arb.)")
    axes[1].set_title("Seebeck vs temperature")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, which='both', alpha=0.3)

    save("fig7_T_dependence.png")
    plt.close()


# =============================================================================
# Figure 8: Transport regimes diagram
# =============================================================================
def figure_8_regimes():
    """Schematic/computational diagram of transport regimes in (W, T) space."""
    W_values = np.logspace(0.5, 4.5, 200)
    T_values = np.linspace(2, 50, 150)

    # We compute ratios:
    #   ell_G / W        -> Gurzhi ratio
    #   ell_5 / W        -> chiral Gurzhi ratio
    # and the ratio of sigma_hydro(W) / sigma_hydro_bulk which indicates
    # how confined the flow is.
    mu = 5.0
    # Fixed tau_ee, scan W and T. tau_mr fixed. tau_5 fixed.
    # ell_G depends on T through w = (4/3) eps, but weakly; we recompute for each T.
    reg = np.zeros((len(T_values), len(W_values)), dtype=float)
    for i, T in enumerate(T_values):
        for j, W in enumerate(W_values):
            p = WeylParams(B=0.0, b_shift=0.0, mu=mu, T=T, W=W)
            sc = hydrodynamic_scales(p)
            r_G = sc['ell_G'] / W
            r_5 = sc['ell_5'] / W
            # reg value encodes regime:
            #   r_G > 1 : viscous Gurzhi matters (Poiseuille)
            #   r_5 > 1 : anomalous Gurzhi matters (chiral Gurzhi)
            # We plot (log10(r_G) + log10(r_5))/2 for color, and overlay
            # the r_G = 1 and r_5 = 1 lines.
            reg[i, j] = np.log10(r_G)

    fig, ax = plt.subplots(1, 1, figsize=(5.5, 3.6))
    X, Y = np.meshgrid(W_values, T_values)
    pcm = ax.pcolormesh(X, Y, reg, cmap='RdBu_r', vmin=-2, vmax=2, shading='auto')
    ax.set_xscale('log')
    ax.set_xlabel(r"$W$ (nm)")
    ax.set_ylabel(r"$T$ (K)")
    ax.set_title("Transport regimes:  color = $\\log_{10}(\\ell_G/W)$")
    cb = plt.colorbar(pcm, ax=ax)
    cb.set_label(r"$\log_{10}(\ell_G / W)$")

    # Overlay the W = 2 ell_G (r_G = 1/2 -> log = -0.3) and W = 2 ell_5 lines
    # Draw for each T
    ellG_line = []; ell5_line = []
    for T in T_values:
        p = WeylParams(B=0.0, b_shift=0.0, mu=mu, T=T, W=100)
        sc = hydrodynamic_scales(p)
        ellG_line.append(2 * sc['ell_G'])
        ell5_line.append(2 * sc['ell_5'])
    ax.plot(ellG_line, T_values, 'k-', lw=2, label=r"$W = 2\ell_G$  (viscous Gurzhi)")
    ax.plot(ell5_line, T_values, 'k--', lw=2, label=r"$W = 2\ell_5$  (chiral Gurzhi)")
    ax.legend(fontsize=9, loc='lower left')

    # Annotate regimes
    ax.annotate("viscous\n(Gurzhi)", xy=(5, 45), fontsize=9, color='white', ha='left')
    ax.annotate("crossover", xy=(2e3, 35), fontsize=9, color='white', ha='center')
    ax.annotate("bulk\nhydro", xy=(2e4, 10), fontsize=9, color='black', ha='center')

    save("fig8_regimes.png")
    plt.close()


# =============================================================================
# Figure 9: Benchmark figure -- shows our result matches known limits
# =============================================================================
def figure_9_benchmark_summary():
    """A single summary figure showing our result agrees with known limits."""
    W_values = np.logspace(0, 4.5, 150)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2))

    # Panel A: viscous Gurzhi vs analytic
    p = WeylParams(B=0.0, b_shift=0.0)
    sigmaH, g_G_ana = [], []
    for W in W_values:
        p.W = W
        out = channel_transport(p)
        sigmaH.append(out['sigma_hydro'])
    sigmaH = np.array(sigmaH)
    sigmaH_bulk = sigmaH[-1]
    sc = hydrodynamic_scales(p)
    g_G_ana = g_gurzhi(W_values / (2 * sc['ell_G']))
    axes[0].semilogx(W_values, sigmaH / sigmaH_bulk, 'b-', lw=2.5, label='numerical')
    axes[0].semilogx(W_values, g_G_ana, 'k--', lw=1.3, label=r'analytic $g(W/2\ell_G)$')
    axes[0].set_xlabel(r"$W$ (nm)")
    axes[0].set_ylabel(r"$\sigma_{\rm hydro}/\sigma_{\rm hydro}^\infty$")
    axes[0].set_title("(a) Gurzhi benchmark")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, which='both', alpha=0.3)

    # Panel B: anomalous Gurzhi vs g_5
    p = WeylParams(B=1.0, b_shift=0.0)
    sigmaA = []
    for W in W_values:
        p.W = W
        out = channel_transport(p)
        sigmaA.append(out['sigma_anom'])
    sigmaA = np.array(sigmaA)
    sigmaA_bulk = sigmaA[-1]
    g5_ana = g_gurzhi(W_values / (2 * hydrodynamic_scales(p)['ell_5']))
    axes[1].semilogx(W_values, sigmaA / sigmaA_bulk, 'r-', lw=2.5, label='numerical')
    axes[1].semilogx(W_values, g5_ana, 'k--', lw=1.3, label=r'analytic $g(W/2\ell_5)$')
    axes[1].set_xlabel(r"$W$ (nm)")
    axes[1].set_ylabel(r"$\sigma_{\rm anom}/\sigma_{\rm anom}^\infty$")
    axes[1].set_title("(b) Sukhachov-Trauzettel benchmark")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, which='both', alpha=0.3)

    # Panel C: negative MR at bulk
    p = WeylParams(W=1e10, b_shift=0.0)
    B_range = np.linspace(0, 1.5, 40)
    sigma_B = []
    for B in B_range:
        p.B = B
        out = channel_transport(p)
        sigma_B.append(out['sigma_xx'])
    sigma_B = np.array(sigma_B)
    sigma_0 = sigma_B[0]
    # Also fit B^2 correction
    ratio = sigma_B / sigma_0
    axes[2].plot(B_range, ratio, 'g-', lw=2.5, label='numerical')
    axes[2].plot(B_range, 1 + (ratio[-1]-1)/(B_range[-1]**2) * B_range**2,
                 'k--', lw=1.3, label=r'$\propto B^2$ fit')
    axes[2].set_xlabel(r"$B$ (dim-less)")
    axes[2].set_ylabel(r"$\sigma_{xx}(B)/\sigma_{xx}(0)$")
    axes[2].set_title("(c) Negative MR (chiral anomaly)")
    axes[2].legend(fontsize=9)
    axes[2].grid(True, which='both', alpha=0.3)

    save("fig9_benchmarks.png")
    plt.close()


if __name__ == "__main__":
    print("Generating figures...")
    figure_1_gurzhi_factor()
    figure_2_viscous_and_chiral_gurzhi()
    figure_3_thermoelectric_gurzhi()
    figure_4_L_contour()
    figure_5_anom_contribution()
    figure_6_poiseuille_profiles()
    figure_7_T_dependence()
    figure_8_regimes()
    figure_9_benchmark_summary()
    print("All figures written.")
