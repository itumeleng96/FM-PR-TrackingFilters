"""Render Log-Likelihood equations — style matches adaptive-filter equation card."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex":       False,
    "mathtext.fontset":  "cm",       # Computer-Modern — matches the reference
    "font.family":       "serif",
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "savefig.facecolor": "white",
})

# Same proportions as the reference card (roughly 3:4 portrait-ish, or wider)
fig, ax = plt.subplots(figsize=(6.5, 5.2))
ax.set_axis_off()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

BLACK = "#111111"
GRAY  = "#444444"

def tx(x, y, s, fs=13, color=BLACK, ha="left", va="center", **kw):
    ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va,
            transform=ax.transAxes, **kw)

X0 = 0.06   # left margin (same feel as reference)

# ── Section 1: Log-Likelihood ───────────────────────────────────────────────
tx(X0, 0.92, "Log-Likelihood", fs=13, fontweight="bold")

tx(X0, 0.80,
   r"$\ell_k \;=\; -\,\dfrac{1}{2}"
   r"\!\left[\;\ln|S_k|\;+\;e_k^{\,T}\,S_k^{-1}\,e_k\;+\;n\ln(2\pi)\;\right]$",
   fs=14.5)

# ── Section 2: Innovation covariance ───────────────────────────────────────
tx(X0, 0.63, "Innovation covariance", fs=13, fontweight="bold")

tx(X0, 0.53,
   r"$S_k \;=\; H_k\,P_k^{-}\,H_k^{T} \;+\; R_k$",
   fs=13.5)

# ── Section 3: Innovation ───────────────────────────────────────────────────
tx(X0, 0.40, "Innovation (residual)", fs=13, fontweight="bold")

tx(X0, 0.30,
   r"$e_k \;=\; z_k \;-\; H_k\,\hat{X}_k^{-}$",
   fs=13.5)

# ── Section 4: Interpretation ───────────────────────────────────────────────
tx(X0, 0.17, "Interpretation", fs=13, fontweight="bold")

tx(X0, 0.08,
   r"Higher (less negative) $\ell_k$  =  better consistency,"
   r"     $P_k^- \downarrow \;\Rightarrow\; \ell_k \to 0$",
   fs=12, color=GRAY)

# ── save ───────────────────────────────────────────────────────────────────
base = ("/home/itumeleng/Documents/Academia/Conference Papers/"
        "CIE_Conference_2026/Presentation/Figures/ll_equations")
fig.savefig(base + ".svg", format="svg", bbox_inches="tight", dpi=300)
fig.savefig(base + ".png", format="png", bbox_inches="tight", dpi=300)
print(f"Saved → {base}.svg / .png")
