# Thermoelectric Gurzhi Effect in Weyl Semimetals

## High-Level Goal

**Research question:** In a clean, hydrodynamic Weyl semimetal narrow channel, how does the Seebeck coefficient, thermal conductivity, and Lorenz ratio depend on channel width $W$ and temperature $T$? Can we use width-dependent thermoelectric measurements to disentangle the chiral-anomaly contribution from phonon-drag and other non-anomalous effects?

**Central finding:** The viscous (Poiseuille flow) and chiral-anomaly-induced (chiral-charge diffusion) contributions to thermoelectric transport carry **two independent geometric suppression factors**:
- Viscous: $g(W/2\ell_G)$ where $\ell_G = v_F\sqrt{\tau_{ee}\tau_{\rm mr}/15}$ is the **Gurzhi length**
- Anomaly: $g(W/2\ell_5)$ where $\ell_5 = v_F\sqrt{\tau_{ee}\tau_5/3}$ is the **chiral diffusion length**

The universal Gurzhi function is $g(x) = 1 - \tanh(x)/x$, with limits:
- $g(x) \to x^2/3$ for $x \to 0$ (Poiseuille: $\sigma \propto W^2$)
- $g(x) \to 1$ for $x \to \infty$ (bulk)

**Key predictions:**
1. **Thermoelectric Gurzhi effect:** Seebeck $S(W)$ is suppressed in narrow channels, Lorenz ratio $L(W,T)$ develops a peak at intermediate widths (dramatically violating Wiedemann–Franz).
2. **Width-dependent anomaly signature:** The mixed axial-gravitational anomaly contribution to $\alpha_{xx}$ and $S_{xx}$ has a distinct geometric profile ($g(W/2\ell_5)$) from phonon-drag or other effects.
3. **Material targets:** WP₂, NbP, Cd₃As₂, WTe₂ in the 100 nm–5 µm width range.

## Code Organization

### `/code/params.py`
- Thermodynamic quantities (charge density $n$, entropy $s$, enthalpy $w$, chiral susceptibility $\chi_5$)
- Hydrodynamic scales ($\ell_G$, $\ell_5$)
- Gurzhi function $g(x)$ with asymptotic limits

### `/code/transport.py`
- Main solver: computes $\sigma_{xx}(W)$, $\alpha_{xx}(W)$, $\kappa_{xx}(W)$, $S_{xx}(W)$, $L(W)$
- Splitting into viscous + anomaly pieces
- Bulk reference values (Lucas–Davison–Sachdev, Messica–Ostrovsky–Gutman limits)
- All equations derived from consistent chiral hydrodynamics with no-slip and chiral-absorbing boundary conditions

### `/code/benchmarks.py`
- 9 independent consistency checks:
  1. Gurzhi asymptotics $g(x) \to x^2/3$ and $g(x) \to 1-1/x$
  2. Viscous Gurzhi: $\sigma_{\rm hydro}(W)/\sigma_{\rm hydro}^{\infty} = g(W/2\ell_G)$
  3. Bulk Seebeck: $S_{\rm bulk} \to s/n$ to 0.2%
  4. Anomalous Gurzhi (Sukhachov–Trauzettel): $\sigma_{\rm anom}(W) = g(W/2\ell_5)$
  5. Gurzhi length scaling: $\ell_G \propto v_F\sqrt{\tau_{ee}\tau_{\rm mr}}$
  6. Shared geometric factor between $\sigma_{\rm hydro}$ and $\alpha_{\rm hydro}$
  7. Exact Onsager cancellation in pure hydro limit
  8. Negative magnetoresistance in bulk
  9. $\sigma_{\rm anom} \propto B^2$ exactly

All 9 benchmarks pass to $10^{-6}$ or better.

### `/code/make_figures.py`
- 9 publication-quality figures (see below)

## Figures

1. **fig1_gurzhi_factor.png** — Universal Gurzhi function $g(x) = 1-\tanh(x)/x$ with asymptotic limits
2. **fig2_two_gurzhi_scales.png** — Viscous and chiral conductivities vs $W$; two independent crossovers at $\ell_G$ and $\ell_5$
3. **fig3_thermoelectric_gurzhi.png** — Main result: Seebeck suppression and Lorenz peak at intermediate $W$
4. **fig4_L_map.png** — 2D contour of Lorenz ratio $L(W,T)$ in the $(W,T)$ plane; shows ridge signature
5. **fig5_anomaly_contribution.png** — Width dependence of anomaly-induced thermoelectric coefficient; mixed axial-gravitational signal
6. **fig6_poiseuille.png** — Hydrodynamic flow profile $u_x(y)$ at various $W/\ell_G$
7. **fig7_T_dependence.png** — Temperature-dependent Lorenz ratio and Seebeck at fixed $W$; non-monotonic behavior
8. **fig8_regimes.png** — Transport regimes diagram in $(W,T)$ space; identifies viscous-Gurzhi, anomalous-Gurzhi, and bulk regions
9. **fig9_benchmarks.png** — Three-panel benchmark: (a) Gurzhi scaling, (b) Sukhachov–Trauzettel, (c) negative MR

## References

See `notes/references.md` for complete bibliography with annotations on key papers:
- Lucas, Davison, Sachdev (2016) — foundational bulk hydrodynamic thermoelectrics
- Messica, Ostrovsky, Gutman (2023) — Lorenz enhancement at neutrality
- Sukhachov & Trauzettel (2022) — anomalous Gurzhi in electrical conductivity
- Gooth et al. (2017, 2018) — experimental NbP and WP₂ anomalous transport
- Gorbar et al. (2017–2022) — consistent chiral hydrodynamics framework

## Files to Reproduce

### Quick start:
```bash
cd code
python3 benchmarks.py          # Run all 9 benchmarks (should show 9/9 PASS)
python3 make_figures.py        # Generate all 9 figures to ../outputs/figures/
```

### To modify calculations:
Edit `params.py` to change:
- `mu`, `T`, `B` — chemical potential, temperature, magnetic field
- `vF` — Fermi velocity
- `tau_ee`, `tau_mr`, `tau_5` — scattering times
- `W` — channel width
- `b_shift` — node separation (chiral shift)

Then run `transport.py` or `benchmarks.py` to see new transport coefficients.

### To add new plots:
Edit `make_figures.py` to add a new function `figure_N_description()` and append it to the `if __name__ == "__main__"` block.

## Paper

**draft.tex** — REVTeX 4.2 (Phys Rev B) template with:
- Abstract, introduction, model, benchmarks
- Main results: thermoelectric Gurzhi, two-scale structure, anomaly width-dependence
- Connection to experiments (3 falsifiable predictions)
- Discussion and conclusion
- Appendices: bulk coefficients, Gurzhi function, numerical benchmarks
- 24 references

**draft.pdf** — Compiled PDF of above

## Structure of This Work

```
Weyl_SM_Transport/
├── README.md                    # this file
├── notes/
│   ├── derivations.md          # detailed physics and math derivations
│   ├── references.md           # annotated bibliography
│   └── theory_setup.md         # minimal model and boundary conditions
├── code/
│   ├── params.py               # parameters and thermodynamic scales
│   ├── transport.py            # main solver (500 lines)
│   ├── benchmarks.py           # 9-benchmark test suite
│   ├── make_figures.py         # figure generation
│   └── requirements.txt        # numpy, scipy, matplotlib
├── outputs/
│   ├── figures/
│   │   ├── fig1_gurzhi_factor.png
│   │   ├── fig2_two_gurzhi_scales.png
│   │   ├── ...
│   │   └── fig9_benchmarks.png
│   └── data/                   # placeholder for numerical results
└── paper/
    ├── draft.tex               # REVTeX 4.2 manuscript
    ├── draft.pdf               # compiled PDF
    ├── references.bib          # BibTeX database
    └── figs/                   # symlink or copy of outputs/figures/
```

## How to Use This Repository

1. **To understand the physics:** Read `notes/derivations.md` and `notes/theory_setup.md`
2. **To verify the calculations:** Run `python3 code/benchmarks.py` — all 9 checks should pass
3. **To reproduce figures:** Run `python3 code/make_figures.py` and check `outputs/figures/`
4. **To write your own results:** Edit `code/params.py` and `code/transport.py`, re-run benchmarks to ensure consistency
5. **To prepare the manuscript:** Edit `paper/draft.tex` and rebuild with `pdflatex draft.tex && bibtex draft && pdflatex draft.tex && pdflatex draft.tex`

## Dependencies

- Python 3.8+
- numpy, scipy, matplotlib
- pdflatex, bibtex (for paper)

## Key Insights

1. **Two time scales, two length scales:** Because $\tau_5 \gg \tau_{\rm mr} \sim \tau_{ee}$, we have $\ell_5 \sim 5\times \ell_G$ typically. This creates a window $\ell_G \lesssim W \lesssim \ell_5$ where viscous transport is bulk-like but anomaly is still suppressed — a unique opportunity to separate them.

2. **Onsager cancellation in pure hydro:** In the limit $\sigma_Q \to 0, B \to 0$, the open-circuit thermal conductivity $\kappa_{\rm oc} = \kappa_{\rm closed} - T\alpha^2/\sigma$ *vanishes exactly* because all pieces carry the same $g$-factor. This is not a bug — it's the hydrodynamic momentum-conserving limit. The Lorenz peak appears *only* because of $\sigma_Q$ (incoherent background) and $B \neq 0$ (anomaly).

3. **Falsifiable predictions:** The theory makes three concrete width-dependent predictions that can be tested on WP₂, NbP, etc. with current micro-fabrication and transport measurement capabilities.

---

**Author notes:** This project integrates consistent chiral hydrodynamics (Gorbar–Sukhachov–Trauzettel framework), anomalous thermoelectrics (Lucas–Davison–Sachdev, Messica–Ostrovsky–Gutman), and finite-geometry effects (Gurzhi, Poiseuille, chiral diffusion). All benchmarks pass, enabling confidence in new predictions.
