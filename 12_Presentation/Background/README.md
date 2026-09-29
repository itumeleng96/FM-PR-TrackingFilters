# Background — FM Illuminator Map

Generates the SA FM broadcast transmitter map showing IcASA-registered FM illuminators of opportunity.

## Pipeline

```
icasa_fm_tables.docx  →  ExtractFM_Sites.py  →  sa_fm_sites.csv
za_provinces.geojson  ┐
sa_fm_sites.csv       →  PlotFM_Illuminators.py  →  Background/Figures/sa_fm_illuminators.{svg,png}
```

## Files

| File | Description |
|---|---|
| `icasa_fm_tables.docx` | Source: ICASA *Optimisation of FM Sound Broadcasting* (2022), Annex tables |
| `za_provinces.geojson` | GADM 4.1 SA province boundaries (level 1) |
| `ExtractFM_Sites.py` | Parses `icasa_fm_tables.docx` → `sa_fm_sites.csv` and `sa_fm_transmitters.csv` |
| `PlotFM_Illuminators.py` | Plots FM sites with ERP service-area rings → `Background/Figures/sa_fm_illuminators.*` |
| `sa_fm_transmitters.csv` | Per-allocation rows (1 759 entries, 87.5–108 MHz) |
| `sa_fm_sites.csv` | Unique sites aggregated by location (597 sites, max ERP per site) |

## Run

```bash
cd 12_Presentation/Background
python3 ExtractFM_Sites.py       # generates sa_fm_sites.csv
python3 PlotFM_Illuminators.py   # writes Figures/sa_fm_illuminators.svg + .png
```

## Data source

ICASA — *Discussion Document on the Optimisation of the FM Sound Broadcasting Band*, April 2022.  
Table 4 (detailed per-station Annex): 1 759 FM allocations across 597 unique transmitter sites.

## Notes

- ERP rings: radius ∝ √ERP (kW), clipped to [8 km, 90 km].
- `za_provinces.geojson` property key tried in order: `"name"` then `"NAME_1"`.
- ATV/DTT broadcast site data (discarded) lives in `Other/`.
