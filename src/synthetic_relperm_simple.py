"""
Simple deterministic synthetic Water-Oil Relative Permeability generator
--------------------------------------------------------------------------
Inputs available : Porosity (phi) and Water Saturation (Sw) from open-hole
                    logs ONLY. No absolute permeability (Ka), no free-water
                    level (FWL).

Calibration data  : Result_Summary_KwKo_YKN.xlsx (4 HFU core samples with
                    Ka, phi, Swir, and full Krw/Kro vs Sw curves)

Method:
 1) Fit a Timur-type permeability transform: log(Ka) = f(log(phi), log(Swir))
    using the 4 core samples (Ka, phi, Swir all known from core).
 2) Fit a Corey model (Krw, Kro vs Sw) to each of the 4 core samples.
 3) Correlate the Corey parameters (Swir, Sor, Krw_max, Kro_max, nw, no)
    against log(Ka) across the 4 samples.
 4) For a NEW sample with only (phi, Sw) known:
      - estimate Ka via the Timur transform (Sw is treated as ~Swir,
        i.e. pick a log depth that looks to be at/near irreducible
        saturation, not in an active transition zone)
      - predict Corey parameters from that Ka via the step-3 correlation
      - build the full synthetic Krw-Kro vs Sw curve

No uncertainty bands - single deterministic curve, as requested.
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# 1. Core calibration data (Ka, phi, Swir, and Krw/Kro vs Sw)
# ---------------------------------------------------------------------
core_data = {
 'HFU-1': dict(Ka=721.799, phi=0.29302, Swir=0.42056,
   Sw=[0.42056,0.486097,0.564864,0.635708,0.703297,0.734614,0.751563,0.772772,0.777936],
   Krw=[0,0.015305,0.024838,0.031163,0.046510,0.053984,0.058545,0.063980,0.064817],
   Kro=[0.344291,0.124977,0.051013,0.020086,0.006302,0.002272,0.001157,0.000152,0]),
 'HFU-2': dict(Ka=198.285, phi=0.28535, Swir=0.43448,
   Sw=[0.43448,0.605054,0.653383,0.719527,0.744471,0.750246,0.759545],
   Krw=[0,0.009076,0.016380,0.028406,0.034004,0.037113,0.037130],
   Kro=[0.241053,0.020088,0.012241,0.001754,0.000287,0.000157,0]),
 'HFU-3': dict(Ka=40.513, phi=0.27884, Swir=0.49115,
   Sw=[0.49115,0.536825,0.592141,0.681924,0.707199,0.786725,0.792319],
   Krw=[0,0.000126,0.000282,0.000895,0.001100,0.001708,0.001755],
   Kro=[0.003044,0.001783,0.001097,0.000339,0.000226,0.0000148,0]),
 'HFU-4': dict(Ka=3.816, phi=0.20806, Swir=0.54385,
   Sw=[0.54385,0.569098,0.610657,0.633551,0.669764,0.703872,0.717129,0.734779],
   Krw=[0,0.001206,0.002588,0.002493,0.004407,0.005253,0.006693,0.007663],
   Kro=[0.057385,0.010852,0.002291,0.001078,0.000626,0.000148,0.0000772,0]),
}

# ---------------------------------------------------------------------
# 2. Timur-type transform: Ka from (phi, Swir) - no permeability needed as input
# ---------------------------------------------------------------------
Ka_arr   = np.array([d['Ka'] for d in core_data.values()])
phi_arr  = np.array([d['phi'] for d in core_data.values()])
Swir_arr = np.array([d['Swir'] for d in core_data.values()])

X = np.column_stack([np.log(phi_arr), np.log(Swir_arr), np.ones(len(Ka_arr))])
y = np.log(Ka_arr)
TIMUR_A, TIMUR_C, TIMUR_D = np.linalg.lstsq(X, y, rcond=None)[0]

def estimate_Ka(phi, Sw):
    """Sw is treated as ~Swir - pick a log point at/near irreducible saturation."""
    return np.exp(TIMUR_A*np.log(phi) + TIMUR_C*np.log(Sw) + TIMUR_D)

# ---------------------------------------------------------------------
# 3. Fit Corey model to each core sample
# ---------------------------------------------------------------------
def fit_corey(d):
    Sw = np.array(d['Sw']); Krw = np.array(d['Krw']); Kro = np.array(d['Kro'])
    Swir, Sor = Sw[0], 1 - Sw[-1]
    Krw_max, Kro_max = Krw[-1], Kro[0]
    Swn = np.clip((Sw - Swir) / (1 - Swir - Sor), 0, 1)
    nw, _ = curve_fit(lambda s, n: Krw_max*s**n, Swn[Krw > 0], Krw[Krw > 0], p0=[2])
    no, _ = curve_fit(lambda s, n: Kro_max*(1-s)**n, Swn[Kro > 0], Kro[Kro > 0], p0=[2])
    return dict(Ka=d['Ka'], Swir=Swir, Sor=Sor, Krw_max=Krw_max, Kro_max=Kro_max,
                nw=nw[0], no=no[0])

fitted = {hfu: fit_corey(d) for hfu, d in core_data.items()}

# ---------------------------------------------------------------------
# 4. Correlate Corey parameters vs log(Ka)
# ---------------------------------------------------------------------
logKa = np.log10([f['Ka'] for f in fitted.values()])
param_names = ['Swir', 'Sor', 'Krw_max', 'Kro_max', 'nw', 'no']
corr_models = {}
for p in param_names:
    vals = np.array([f[p] for f in fitted.values()])
    corr_models[p] = np.polyfit(logKa, vals, 1)   # (slope, intercept)

def predict_corey_params(Ka):
    lk = np.log10(Ka)
    return {p: a*lk + b for p, (a, b) in corr_models.items()}

# ---------------------------------------------------------------------
# 5. MAIN FUNCTION: phi, Sw (log) -> synthetic Krw-Kro curve (deterministic)
# ---------------------------------------------------------------------
def generate_synthetic_relperm(phi, Sw, n_points=25):
    """
    phi : porosity from log (fraction)
    Sw  : water saturation from log (fraction) - should be picked from a
          depth judged to be at/near irreducible saturation.

    Returns: Sw_grid, Krw, Kro, Ka_estimated, corey_params (dict)
    """
    Ka_est = estimate_Ka(phi, Sw)
    p = predict_corey_params(Ka_est)

    Swir, Sor = p['Swir'], p['Sor']
    Sw_grid = np.linspace(Swir, 1 - Sor, n_points)
    Swn = np.clip((Sw_grid - Swir) / (1 - Swir - Sor), 0, 1)
    Krw = p['Krw_max'] * Swn**p['nw']
    Kro = p['Kro_max'] * (1 - Swn)**p['no']
    return Sw_grid, Krw, Kro, Ka_est, p

# ---------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------
if __name__ == "__main__":
    phi_input, Sw_input = 0.25, 0.48   # <-- replace with your log values

    Sw_grid, Krw, Kro, Ka_est, params = generate_synthetic_relperm(phi_input, Sw_input)

    print(f"Input: phi={phi_input}, Sw={Sw_input}")
    print(f"Estimated Ka: {Ka_est:.2f} mD")
    print("Corey parameters:")
    for k, v in params.items():
        print(f"  {k}: {v:.4f}")
    print("\nSw       Krw        Kro")
    for s, w, o in zip(Sw_grid, Krw, Kro):
        print(f"{s:.4f}   {w:.6f}   {o:.6f}")

    fig, ax = plt.subplots(figsize=(8,6))

    # --- underlying core data points used to build the correlations ---
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(core_data)))
    for (hfu, d), c in zip(core_data.items(), colors):
        ax.plot(d['Sw'], d['Krw'], 'o', color=c, ms=6, mec='k', mew=0.5,
                label=f"{hfu} core data (Ka={d['Ka']:.0f}mD)")
        ax.plot(d['Sw'], d['Kro'], 's', color=c, ms=6, mec='k', mew=0.5)

        # light fitted Corey curve for that core sample, for visual continuity
        f = fitted[hfu]
        Swn_line = np.linspace(0, 1, 60)
        Sw_line = f['Swir'] + Swn_line*(1 - f['Swir'] - f['Sor'])
        ax.plot(Sw_line, f['Krw_max']*Swn_line**f['nw'], '-', color=c, lw=0.8, alpha=0.5)
        ax.plot(Sw_line, f['Kro_max']*(1-Swn_line)**f['no'], '--', color=c, lw=0.8, alpha=0.5)

    # --- the new synthetic curve, highlighted ---
    ax.plot(Sw_grid, Krw, 'b-o', ms=5, lw=2.5, label=f'Krw synthetic (phi={phi_input}, Sw={Sw_input})', zorder=5)
    ax.plot(Sw_grid, Kro, 'r-s', ms=5, lw=2.5, label=f'Kro synthetic (Ka est.={Ka_est:.0f}mD)', zorder=5)

    ax.set_xlabel('Sw (fraction)'); ax.set_ylabel('Kr (fraction)')
    ax.set_title('Synthetic curve vs. the 4 core samples used to build the correlations\n(o = Krw, s = Kro; dashed/thin lines = fitted Corey curves for core data)')
    ax.legend(fontsize=7, loc='upper center', ncol=2)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('/home/rian/python_project/myvenv/reservoir_analysis/images/synthetic_relperm_simple.png', dpi=150)
    print("\nSaved plot to synthetic_relperm_simple.png")
