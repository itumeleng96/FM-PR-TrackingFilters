"""
Extract unique broadcast transmitter sites from ICASA Terrestrial Broadcasting
Frequency Plan (icasa_tbfp.xlsx). Combines Annex E (ATV) and Annex G (DTT),
converts DMS lat/long to decimal, and writes a deduped CSV.

Author: Itumeleng Malemela
"""

import csv
import re
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent
XLSX = HERE / "icasa_tbfp.xlsx"
OUT_CSV = HERE / "sa_broadcast_sites.csv"


DMS_RE = re.compile(r"^\s*(\d+)\s*([NSEW])\s*(\d+)\s*(\d+)?\s*$", re.IGNORECASE)


def dms_to_decimal(s):
    """Parse ICASA lat/long like '33S14 00' or '19E17 14' to decimal degrees."""
    if s is None:
        return None
    m = DMS_RE.match(str(s))
    if not m:
        return None
    deg = int(m.group(1))
    hemi = m.group(2).upper()
    minutes = int(m.group(3))
    seconds = int(m.group(4)) if m.group(4) else 0
    dec = deg + minutes / 60.0 + seconds / 3600.0
    if hemi in ("S", "W"):
        dec = -dec
    return round(dec, 5)


def parse_sheet(ws, service_tag, freq_col_idx, erp_col_idx):
    """Yield (station, lat, lon, freq_mhz, erp_kw, service) per data row."""
    header_row = None
    for i, row in enumerate(ws.iter_rows(values_only=True), start=1):
        # ICASA rows have first cell = integer NO for data rows
        if isinstance(row[0], int):
            station = row[1]
            lat_dms = row[2]
            lon_dms = row[3]
            freq = row[freq_col_idx]
            erp = row[erp_col_idx]
            if station is None or lat_dms is None or lon_dms is None:
                continue
            lat = dms_to_decimal(lat_dms)
            lon = dms_to_decimal(lon_dms)
            if lat is None or lon is None:
                continue
            yield (str(station).strip(), lat, lon, freq, erp, service_tag)


def main():
    wb = openpyxl.load_workbook(XLSX, data_only=True)

    rows = []
    # Column indices are 0-based, matching the header we saw:
    # ATV:  NO | STATION | LAT | LONG | FREQ | CH | OFFSET | ERP | POL | ...
    # DTT:  NO | STATION | LAT | LONG | FREQ | CH | ERP    | POL | ...
    rows += list(parse_sheet(wb["AnnexE ATV"], "ATV", freq_col_idx=4, erp_col_idx=7))
    rows += list(parse_sheet(wb["AnnexG DTT"], "DTT", freq_col_idx=4, erp_col_idx=6))

    print(f"raw allocation rows: {len(rows)}")

    # Aggregate to unique sites
    sites = {}
    for station, lat, lon, freq, erp, service in rows:
        key = (station, round(lat, 3), round(lon, 3))
        s = sites.setdefault(
            key,
            {
                "station": station,
                "lat": lat,
                "lon": lon,
                "atv_count": 0,
                "dtt_count": 0,
                "max_erp_kw": 0.0,
            },
        )
        if service == "ATV":
            s["atv_count"] += 1
        else:
            s["dtt_count"] += 1
        try:
            erp_val = float(erp)
            if erp_val > s["max_erp_kw"]:
                s["max_erp_kw"] = erp_val
        except (TypeError, ValueError):
            pass

    unique = sorted(sites.values(), key=lambda x: x["station"])
    print(f"unique sites: {len(unique)}")

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["station", "lat", "lon", "atv_count", "dtt_count", "max_erp_kw"]
        )
        w.writeheader()
        w.writerows(unique)
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
