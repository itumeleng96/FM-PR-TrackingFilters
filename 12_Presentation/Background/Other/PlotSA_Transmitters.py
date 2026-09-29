"""
Plot ICASA-registered SA broadcast transmitter sites over provincial map.
Outputs Figures/sa_broadcast_sites.svg + .png.

Requires:
  sa_broadcast_sites.csv   (from extract_sites.py)
  gadm41_ZAF_1.json        (GADM 4.1 South Africa provinces GeoJSON)

Author: Itumeleng Malemela
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "sa_broadcast_sites.csv"
GEOJSON = HERE / "gadm41_ZAF_1.json"
OUT_DIR = HERE.parent / "Figures"
OUT_DIR.mkdir(exist_ok=True)

# Prieska Tx from FERS scenario_1_singleFile.fersxml
# The scenario places the Tx at UTM-local metres; the real Prieska SABC site is
# at approx -29.669, 22.749. We label it as the FERS scenario transmitter.
PRIESKA = {"name": "Prieska Tx (94 MHz)", "lat": -29.669, "lon": 22.749}


def load_sites():
    sites = []
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            row["lat"] = float(row["lat"])
            row["lon"] = float(row["lon"])
            row["atv_count"] = int(row["atv_count"])
            row["dtt_count"] = int(row["dtt_count"])
            row["max_erp_kw"] = float(row["max_erp_kw"] or 0)
            sites.append(row)
    return sites


def load_provinces():
    with open(GEOJSON) as f:
        gj = json.load(f)
    provs = []
    for feat in gj["features"]:
        name = feat["properties"].get("NAME_1", "?")
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

    fig, ax = plt.subplots(figsize=(10, 8), dpi=140)
    fig.patch.set_facecolor("white")

    draw_provinces(ax, provs)

    # Base scatter — all sites
    lons = [s["lon"] for s in sites]
    lats = [s["lat"] for s in sites]
    services = [s["atv_count"] + s["dtt_count"] for s in sites]
    # Size proportional to number of services hosted (clipped)
    sizes = [8 + min(sv, 12) * 3 for sv in services]

    sc = ax.scatter(
        lons, lats, s=sizes,
        c="#1f6feb", alpha=0.55,
        edgecolor="white", linewidth=0.4,
        zorder=3, label="ICASA broadcast sites",
    )

    # Highlight Prieska
    ax.scatter(
        [PRIESKA["lon"]], [PRIESKA["lat"]],
        s=260, marker="*",
        c="#d62728", edgecolor="white", linewidth=1.4,
        zorder=5, label=PRIESKA["name"],
    )
    ax.annotate(
        PRIESKA["name"],
        xy=(PRIESKA["lon"], PRIESKA["lat"]),
        xytext=(PRIESKA["lon"] + 1.0, PRIESKA["lat"] - 1.8),
        fontsize=10, fontweight="bold", color="#8a1c1f",
        arrowprops=dict(arrowstyle="-", color="#8a1c1f", lw=0.9),
        zorder=6,
    )

    # Cosmetic axes
    ax.set_xlim(15.5, 33.5)
    ax.set_ylim(-35.5, -21.5)
    ax.set_aspect(1 / 0.85)  # approx aspect at SA latitudes
    ax.set_xlabel("Longitude [°E]", fontsize=11)
    ax.set_ylabel("Latitude [°S]", fontsize=11)
    ax.set_title(
        "ICASA-registered terrestrial broadcast transmitter sites\n"
        "(FM 87.5–108 MHz co-located on same Sentech infrastructure)",
        fontsize=12, fontweight="bold",
    )

    # Legend combining size cue
    handles = [
        Line2D([], [], marker="o", color="w",
               markerfacecolor="#1f6feb", markeredgecolor="white",
               markersize=8, label=f"Broadcast site (n = {len(sites)})"),
        Line2D([], [], marker="*", color="w",
               markerfacecolor="#d62728", markeredgecolor="white",
               markersize=15, label=PRIESKA["name"]),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=True,
              facecolor="white", edgecolor="#cccccc", fontsize=10)

    ax.grid(True, color="#dddddd", linewidth=0.5, alpha=0.6, zorder=0)
    for spine in ax.spines.values():
        spine.set_color("#666666")
        spine.set_linewidth(0.8)

    fig.tight_layout()

    svg = OUT_DIR / "sa_broadcast_sites.svg"
    png = OUT_DIR / "sa_broadcast_sites.png"
    fig.savefig(svg, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(png, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"wrote {svg}")
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
