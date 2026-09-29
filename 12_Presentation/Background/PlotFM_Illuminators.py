"""
Plot ICASA-registered FM broadcast transmitter sites over SA provincial map.
Draws indicative service-area rings (radius ∝ √ERP) per site.
Outputs Background/Figures/sa_fm_illuminators.svg + .png.

Requires (same directory):
  sa_fm_sites.csv       -- from ExtractFM_Sites.py
  za_provinces.geojson  -- GADM 4.1 SA provinces GeoJSON (gadm41_ZAF_1.json)

Author: Itumeleng Malemela
"""

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "sa_fm_sites.csv"
GEOJSON = HERE / "za_provinces.geojson"
OUT_DIR = HERE / "Figures"
OUT_DIR.mkdir(exist_ok=True)


def load_sites():
    sites = []
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            row["lat"] = float(row["lat"])
            row["lon"] = float(row["lon"])
            row["fm_count"] = int(row["fm_count"])
            row["max_erp_kw"] = float(row["max_erp_kw"] or 0)
            sites.append(row)
    return sites


def load_provinces():
    with open(GEOJSON) as f:
        gj = json.load(f)
    provs = []
    for feat in gj["features"]:
        name = feat["properties"].get("name", feat["properties"].get("NAME_1", "?"))
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            polys = [geom["coordinates"]]
        elif geom["type"] == "MultiPolygon":
            polys = geom["coordinates"]
        else:
            continue
        provs.append((name, polys))
    return provs


def draw_provinces(ax, provs):
    fill = "#f4f2ec"
    edge = "#8a8577"
    for name, polys in provs:
        for poly in polys:
            outer = poly[0]
            xs = [pt[0] for pt in outer]
            ys = [pt[1] for pt in outer]
            ax.fill(xs, ys, facecolor=fill, edgecolor=edge, linewidth=0.6, zorder=1)


def main():
    sites = load_sites()
    provs = load_provinces()
    print(f"sites: {len(sites)}, provinces: {len(provs)}")

    fig, ax = plt.subplots(figsize=(8.5, 9.5), dpi=140)
    fig.patch.set_facecolor("white")

    draw_provinces(ax, provs)

    # Indicative service-area rings, radius ∝ √ERP (kW).
    # r_km ≈ 12·√ERP, clipped to [8 km, 90 km] so ultra-low-power
    # community sites remain visible and 100 kW giants don't dominate.
    def erp_to_radius_km(erp_kw):
        r = 12.0 * math.sqrt(max(erp_kw, 0.01))
        return max(8.0, min(90.0, r))

    for s in sites:
        r_km = erp_to_radius_km(s["max_erp_kw"])
        r_lon = r_km / (111.0 * math.cos(math.radians(s["lat"])))
        ax.add_patch(Circle(
            (s["lon"], s["lat"]), radius=r_lon,
            facecolor="#1f6feb", edgecolor="none",
            alpha=0.18, zorder=2,
        ))

    lons = [s["lon"] for s in sites]
    lats = [s["lat"] for s in sites]
    ax.scatter(
        lons, lats, s=6,
        c="#0b3b8a", alpha=0.9,
        edgecolor="none", zorder=3,
    )

    # Clean map — no ticks, no spines, no grid
    ax.set_xlim(15.5, 33.5)
    ax.set_ylim(-35.5, -21.5)
    ax.set_aspect(1.0 / math.cos(math.radians(30)))  # equal-km at SA mid-latitudes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Legend: site marker + ERP ring-size key (2×2 grid below map)
    ring_face = (0.121, 0.435, 0.921, 0.35)  # #1f6feb @ alpha 0.35
    site_dot = Line2D(
        [], [], marker="o", color="w",
        markerfacecolor="#0b3b8a", markeredgecolor="none",
        markersize=6, label=f"Broadcast site (n = {len(sites)})",
    )
    r_low = Line2D([], [], marker="o", color="w",
                   markerfacecolor=ring_face, markeredgecolor="none",
                   markersize=10, label="ERP 0.1 kW  (~8 km)")
    r_mid = Line2D([], [], marker="o", color="w",
                   markerfacecolor=ring_face, markeredgecolor="none",
                   markersize=18, label="ERP 10 kW  (~38 km)")
    r_high = Line2D([], [], marker="o", color="w",
                    markerfacecolor=ring_face, markeredgecolor="none",
                    markersize=28, label="ERP 100 kW (~90 km)")
    legend = ax.legend(
        handles=[site_dot, r_low, r_mid, r_high],
        title="Indicative service area (radius ∝ √ERP)",
        loc="upper center", bbox_to_anchor=(0.5, -0.02),
        ncol=2, frameon=True, facecolor="white", edgecolor="#cccccc",
        fontsize=10, title_fontsize=10,
        columnspacing=3.0, handletextpad=1.4,
        labelspacing=1.4, borderpad=1.0,
    )
    legend.get_title().set_fontweight("bold")

    ax.set_title("FM Radio Transmitters in SA", fontsize=13, fontweight="bold", pad=2)

    fig.tight_layout(rect=[0, 0.13, 1, 1])

    svg = OUT_DIR / "sa_fm_illuminators.svg"
    png = OUT_DIR / "sa_fm_illuminators.png"
    fig.savefig(svg, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(png, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"wrote {svg}")
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
