# Other — ATV/DTT Broadcast Sites (Discarded)

Extracts and plots all IcASA-registered terrestrial broadcast sites (ATV + DTT).
Superseded by the FM-only pipeline in the parent directory.

## Files

| File | Description |
|---|---|
| `icasa_tbfp.xlsx` | Source: ICASA *Terrestrial Broadcasting Frequency Plan* — Annex E (ATV) + Annex G (DTT) |
| `ExtractSA_Sites.py` | Parses `icasa_tbfp.xlsx` → `sa_broadcast_sites.csv` (313 unique ATV/DTT sites) |
| `PlotSA_Transmitters.py` | Maps all broadcast sites + Prieska Tx → `Figures/sa_broadcast_sites.*` |

## Notes

- Superseded because FM allocations come from a separate ICASA document (`icasa_fm_tables.docx`).
- FM 87.5–108 MHz shares the same Sentech infrastructure as ATV/DTT, but the FM-specific pipeline (parent `Background/`) uses precise per-frequency ERP values from the FM optimisation tables.
