"""
Extract FM 87.5–108 MHz per-transmitter allocations from ICASA's
'Discussion Document on the Optimisation of the FM Sound Broadcasting - Tables'
(icasa_fm_tables.docx). Uses Table 4 of the document, which is the detailed
per-station Annex.

Output: sa_fm_transmitters.csv with one row per (site, frequency) allocation.
Also emits sa_fm_sites.csv with unique sites (max ERP + count of allocations).

Author: Itumeleng Malemela
"""

import csv
import re
from pathlib import Path

from docx import Document

HERE = Path(__file__).resolve().parent
DOCX = HERE / "icasa_fm_tables.docx"
OUT_ALL = HERE / "sa_fm_transmitters.csv"
OUT_SITES = HERE / "sa_fm_sites.csv"


DMS_RE = re.compile(r"^\s*(\d+)\s*([NSEW])\s*(\d+)\s*(\d+(?:\.\d+)?)?\s*$", re.IGNORECASE)


def dms_to_decimal(s):
    if s is None:
        return None
    m = DMS_RE.match(str(s).strip())
    if not m:
        return None
    deg = int(m.group(1))
    hemi = m.group(2).upper()
    minutes = int(m.group(3))
    seconds = float(m.group(4)) if m.group(4) else 0
    dec = deg + minutes / 60.0 + seconds / 3600.0
    if hemi in ("S", "W"):
        dec = -dec
    return round(dec, 5)


def cell_text(cell):
    return " ".join(cell.text.split()).strip()


def parse_annex_table(table):
    """Yield dicts for each valid data row of Table 4 (the FM Annex)."""
    rows = list(table.rows)
    # Skip the two header rows.
    for row in rows[2:]:
        cells = [cell_text(c) for c in row.cells]
        # Expected 12 columns: NO, STATION, LAT, LONG, FREQ, ERP, POL, PROG, STAT, CAT, ONAIR, COMMENTS
        if len(cells) < 10:
            continue
        no, station, lat_dms, lon_dms, freq, erp, pol, prog, stat, cat = cells[:10]
        if not no.isdigit():
            continue
        lat = dms_to_decimal(lat_dms)
        lon = dms_to_decimal(lon_dms)
        if lat is None or lon is None:
            continue
        try:
            freq_mhz = float(freq)
        except ValueError:
            continue
        try:
            erp_kw = float(erp)
        except ValueError:
            erp_kw = None
        yield {
            "station": station,
            "lat": lat,
            "lon": lon,
            "freq_mhz": freq_mhz,
            "erp_kw": erp_kw,
            "polarisation": pol,
            "program": prog,
            "status": stat,
            "category": cat,
        }


def main():
    doc = Document(DOCX)
    annex = doc.tables[4]  # confirmed 1765 rows, 12 cols
    rows = list(parse_annex_table(annex))
    print(f"parsed FM allocations: {len(rows)}")

    # Filter FM broadcast band 87.5–108 MHz
    rows = [r for r in rows if 87.4 <= r["freq_mhz"] <= 108.1]
    print(f"within FM band (87.5–108 MHz): {len(rows)}")

    with open(OUT_ALL, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT_ALL}")

    # Aggregate to unique sites (station + rounded coords)
    sites = {}
    for r in rows:
        key = (r["station"], round(r["lat"], 3), round(r["lon"], 3))
        s = sites.setdefault(
            key,
            {
                "station": r["station"],
                "lat": r["lat"],
                "lon": r["lon"],
                "fm_count": 0,
                "max_erp_kw": 0.0,
                "programs": set(),
            },
        )
        s["fm_count"] += 1
        if r["erp_kw"] is not None and r["erp_kw"] > s["max_erp_kw"]:
            s["max_erp_kw"] = r["erp_kw"]
        if r["program"]:
            s["programs"].add(r["program"])

    unique = sorted(sites.values(), key=lambda x: x["station"])
    print(f"unique FM sites: {len(unique)}")

    with open(OUT_SITES, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "station", "lat", "lon", "fm_count", "max_erp_kw", "programs"
        ])
        w.writeheader()
        for s in unique:
            w.writerow({
                "station": s["station"],
                "lat": s["lat"],
                "lon": s["lon"],
                "fm_count": s["fm_count"],
                "max_erp_kw": s["max_erp_kw"],
                "programs": "; ".join(sorted(s["programs"])),
            })
    print(f"wrote {OUT_SITES}")


if __name__ == "__main__":
    main()
