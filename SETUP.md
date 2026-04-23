# Weyl Semimetal Transport Project - Setup Instructions

## Quick Start

All files are now ready in `/mnt/user-data/outputs/`. You can download them directly or follow the instructions below to set up your project directory.

## Project Structure

Create this directory structure on your local machine:

```
Weyl_SM_Transport/
├── README.md
├── notes/
│   ├── derivations.md
│   ├── references.md
│   └── theory_setup.md
├── code/
│   ├── params.py
│   ├── transport.py
│   ├── benchmarks.py
│   ├── make_figures.py
│   └── requirements.txt
├── outputs/
│   ├── figures/
│   │   ├── fig1_gurzhi_factor.png
│   │   ├── fig2_two_gurzhi_scales.png
│   │   ├── ... (8 more figures)
│   │   └── fig9_benchmarks.png
│   └── data/
└── paper/
    ├── draft.tex
    ├── draft.pdf
    └── references.bib
```

## Files to Download

### Core Documentation (3 files)
- **README.md** - High-level project overview
- **notes/derivations.md** - Complete physics derivations
- **notes/references.md** - Annotated bibliography (download as references.md)

### Python Code (4 files)
- **code/params.py** - Parameters and thermodynamic scales
- **code/transport.py** - Main hydrodynamic solver
- **code/benchmarks.py** - 9-benchmark validation suite
- **code/make_figures.py** - Figure generation script
- **code/requirements.txt** - Python dependencies

### Figures (9 PNG files)
All images for the paper, ready to embed in LaTeX:
- fig1_gurzhi_factor.png
- fig2_two_gurzhi_scales.png
- fig3_thermoelectric_gurzhi.png
- fig4_L_map.png
- fig5_anomaly_contribution.png
- fig6_poiseuille.png
- fig7_T_dependence.png
- fig8_regimes.png
- fig9_benchmarks.png

### Paper (3 files)
- **paper/draft.tex** - REVTeX manuscript (standard article class)
- **paper/draft.pdf** - Compiled PDF (5 pages)
- **paper/references.bib** - BibTeX reference database

## How to Use

### 1. Setup and Verify
```bash
cd Weyl_SM_Transport/code
pip install -r requirements.txt
python3 benchmarks.py
```

Expected output: `9/9 benchmarks passed`

### 2. Generate Fresh Figures
```bash
python3 make_figures.py
# Outputs to ../outputs/figures/
```

### 3. Customize Calculations

Edit `params.py`:
```python
WeylParams(
    mu=5.0,           # Change chemical potential (meV)
    T=10.0,           # Temperature (K)
    W=500.0,          # Channel width (nm)
    vF=300.0,         # Fermi velocity (nm/ps)
    tau_ee=1.0,       # e-e scattering time (ps)
    tau_mr=10.0,      # momentum relaxing time (ps)
    tau_5=50.0,       # intervalley time (ps)
    B=0.3             # Magnetic field (dimensionless)
)
```

Then re-run benchmarks or generate new figures.

### 4. Edit and Compile Paper
```bash
cd paper/
pdflatex draft.tex
# Edit draft.tex with your text editor
# Recompile:
pdflatex draft.tex
```

For bibliography:
```bash
pdflatex draft.tex
bibtex draft
pdflatex draft.tex
pdflatex draft.tex
```

## File Contents Summary

### README.md (8.7 KB)
- High-level research question and central findings
- Code organization and descriptions
- Key insights and how to use the repository

### derivations.md (20 KB)
- Complete theoretical framework in 10 sections:
  1. Thermodynamic foundation
  2. Hydrodynamic equations in slab geometry
  3. Valley-odd (chiral) charge and anomaly
  4. Linear response and transport coefficients
  5. Thermoelectric coefficients
  6. Two independent geometric scales
  7. Limiting cases and benchmarks
  8. Pure hydrodynamic limit and Onsager cancellation
  9. Connection to experiments
  10. Numerical implementation

### references.md (when needed)
- Annotated bibliography with key papers

### Code Files

**params.py** (7.5 KB):
- WeylParams dataclass
- thermodynamic_quantities() function
- hydrodynamic_scales() function
- g_gurzhi() Gurzhi correction factor
- Self-test on import

**transport.py** (17 KB):
- _coupling_coefficients() - derived couplings
- _drude_like_prefactors() - force-to-velocity conversion
- _chiral_suppression() - chiral-charge geometric factors
- bulk_conductivity() - bulk reference values
- channel_transport() - main solver
- seebeck_bulk_reference() - normalization
- Self-consistency checks

**benchmarks.py** (7.8 KB):
- 9 independent consistency checks:
  1. Gurzhi asymptotics
  2. Viscous Gurzhi scaling
  3. Bulk Seebeck formula
  4. Anomalous Gurzhi
  5. Hydrodynamic scale consistency
  6. Shared geometric factors
  7. Wiedemann-Franz in pure hydro limit
  8. Negative magnetoresistance
  9. B² scaling of anomalous conductivity

**make_figures.py** (19 KB):
- 9 figure generation functions
- All publication-quality plots
- Configurable parameters and styles

### Paper (draft.tex)

5-page manuscript in standard LaTeX article format:
- Abstract
- Introduction (3 paragraphs)
- Model (3 subsections)
- Benchmarks (4 bullet points)
- Results (3 subsections)
- Experimental predictions (3 predictions)
- Discussion
- Conclusion
- References (30 citations)

## Dependencies

### Python
```
numpy>=1.20      # Numerical arrays
scipy>=1.7       # Special functions
matplotlib>=3.5  # Plotting
```

### LaTeX
- pdflatex
- Standard packages (amsmath, amssymb, graphicx, hyperref, xcolor)

## Typical Workflow

```
1. Read README.md to understand the project
2. Read notes/derivations.md for physics background
3. Run code/benchmarks.py to verify setup
4. Edit code/params.py to change parameters
5. Run code/make_figures.py to generate new plots
6. Edit paper/draft.tex with your results
7. Compile with pdflatex to create draft.pdf
```

## Key Physics Points to Remember

1. **Two Independent Gurzhi Scales:**
   - Viscous: ℓ_G = v_F √(τ_ee τ_mr / 15) ≈ 245 nm
   - Chiral: ℓ_5 = v_F √(τ_ee τ_5 / 3) ≈ 1225 nm
   - Ratio: ℓ_5/ℓ_G = √(τ_5/τ_mr) ~ 5

2. **Universal Gurzhi Function:**
   - g(x) = 1 - tanh(x)/x where x = W / (2ℓ)
   - Small x (Poiseuille): g(x) ~ x²/3 → σ ∝ W²
   - Large x (bulk): g(x) → 1

3. **Thermoelectric Gurzhi Effect:**
   - Seebeck suppressed at small W
   - Lorenz number peaks at intermediate W ~ ℓ_G
   - Peak reaches L ~ 3-10 × L₀

4. **Falsifiable Predictions:**
   - (P1) Width-dependent Seebeck in WP₂, NbP
   - (P2) Lorenz ridge in (W,T) plane
   - (P3) Anomaly window where viscous ≈ bulk but anomaly ~ g(W/2ℓ_5)

## Support & Troubleshooting

If benchmarks fail:
```bash
python3 -c "from code.params import WeylParams; from code.transport import channel_transport; p = WeylParams(); out = channel_transport(p); print('Seebeck:', out['S_xx'])"
```

If figures don't generate:
- Check that all PNG output paths exist
- Verify matplotlib is installed: `python3 -c "import matplotlib; print(matplotlib.__version__)"`

If LaTeX won't compile:
- Ensure pdflatex is in PATH: `which pdflatex`
- Try: `pdflatex --version`

---

**All files are now in `/mnt/user-data/outputs/` ready for download and local use.**
