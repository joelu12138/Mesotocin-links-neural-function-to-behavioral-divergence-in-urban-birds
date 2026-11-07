# Mesotocin links neural function to behavioral divergence in urban birds

R code for the statistical and transcriptomic analyses of urban–rural divergence in
Eurasian tree sparrows (*Passer montanus*). Scripts are ordered to follow the
Materials and Methods sections of the manuscript.

## Repository layout (methods order)

| Dir | Methods section | Scripts |
|---|---|---|
| `R/01_behavior` | Behavioral assays | `boldness_ttest`, `spatial_foraging_rm_anova`, `sex_effects_rm_anova`, `social_analysis`, `flocking_analysis`, `flocking_frequency_analysis` |
| `R/02_morphology_histology` | Morphological & histological measurements | `morphology_ttest`, `nissl_analysis` |
| `R/03_spines_synapse` | Dendritic spines & synaptic ultrastructure | `spine_analysis`, `synapse_analysis`, `active_zone_analysis`, `psd_analysis` |
| `R/04_snrna_seq` | Single-nucleus RNA-seq; DE / enrichment / gene-set scoring | `single_cell_diff`, `neuronal_subtype_diff`, `kegg_enrichment` |
| `R/05_network_communication` | Co-expression network; cell–cell communication | `hdwgcna`, `cellchat`, `build_sparrow_cellchat_db` (ships `CellChatDB.treeSparrow_conservative.rda`) |
| `R/06_spatial` | Stereo-seq bin100 regional DE; cellbin RCTD annotation; neuronal regional DE; density; composition chi-square+FDR; decomposition | `spatial_diff_genes`, `spatial_neuron_diff_genes`, `spatial_rctd`, `neuron_subtype_density`, `spatial_composition`, `decomposition_region_subtype` |
| `R/07_neurogenesis_MT` | Adult-born neuronal labeling; regional MT immunoreactivity | `neurogenesis_analysis`, `neurogenesis_1d_analysis`, `mt_ir_analysis`, `mt_ir_sex_analysis` |
| `R/08_electrophysiology` | Electrophysiological recording | `discharge_analysis`, `spontaneous_firing_rm_anova`, `mepsc_rm_anova`, `plot_firing_rate`, `plot_mepsc`, `figs11_descriptive_stats`, `figs11_ttest` |
| `R/09_pharmacology` | In vivo pharmacological interventions | `spatial_cognition_rm_anova`, `social_verification`, `flocking_verification`, `mt_ir_verification_urban`, `mt_ir_verification_rural` (`ota_spatial_cognition` — superseded) |
| `R/10_supplementary` | Supplementary tables | `generate_supp_tables` |

> **Note on `R/06_spatial/decomposition_state_vs_composition_OTHER_PROJECT.R`:**
> this script uses song-system regions (CMM/HVC/LAMN/NCM/RA) and HA/LA groups from a
> *separate* study (cross-scale vocal–auditory mismatch). It is retained here for
> reference only and is **not** part of the analyses reported in this manuscript.
> The decomposition used in this paper is `decomposition_region_subtype.R`
> (HP/NCL/LS/PVN regions, A1/A2 groups).

> **Note on `R/09_pharmacology/ota_spatial_cognition.R`:** superseded. It uses the
> older `Relative_exploration_time` readout; the current analysis uses
> `total_exploration_time` in `spatial_cognition_rm_anova.R`. Retained for
> provenance only.

> **Note on the CellChat database:** `cellchat.R` uses the species-adapted
> `CellChatDB.treeSparrow_conservative` (1-to-1 ortholog map, **1,757**
> ligand-receptor pairs, matching the Methods), shipped as
> `R/05_network_communication/CellChatDB.treeSparrow_conservative.rda` and built
> from `CellChatDB.human` via OrthoFinder by `build_sparrow_cellchat_db.R`. That
> script also produces a full ortholog variant (`CellChatDB.treeSparrow`,
> **12,763** pairs). The build script's ortholog TSV inputs are upstream
> artifacts not included here; the prebuilt `.rda` is ready to use. The DB and
> the expression matrix share the same gene-symbol namespace.

## Dependencies

- R ≥ 4.5
- Statistics: `car`, `emmeans`, `afex`, `lme4`, `lmerTest`, `dplyr`, `tidyr`, `openxlsx`, `ggplot2`
- Transcriptomics: `Seurat` (v5), `harmony`, `irGSEA`, `AUCell`, `UCell`, `clusterProfiler`,
  `hdWGCNA`, `WGCNA`, `CellChat`, `spacexr`, `STdeconvolve`, `scRNAtoolVis`

## Data directory

Input files are **not** versioned (large / external). All scripts read inputs
from and write outputs to a local `data/` directory at the repository root
(relative paths). To run:

1. Place the input CSV / Rdata files in `data/`, using the file names each script
   expects (e.g. `AA.csv`, `Nissl bodies.csv`, `sce_all_integrated.Rdata`,
   `{region}_{group}.bin100.Rdata`, `{region}_{group}.sub.Rdata`).
2. Run scripts from the repository root so `data/` resolves correctly.

Sequencing raw data are deposited separately (see the manuscript Data
Availability statement); `data/` holds the processed tables/objects the scripts
consume.

The gene sets used for irGSEA/AUCell scoring (`marker.csv`, read by
`R/04_snrna_seq/single_cell_diff.R` and `neuronal_subtype_diff.R`) are listed in
Supplementary Table S48 of the manuscript and are not versioned here; save them
as `data/marker.csv` with columns `names` (gene-set name) and `SYMBOL`
(human-ortholog gene symbol) before running those scripts.

## Running

```bash
Rscript R/01_behavior/boldness_ttest.R
```

Output tables are written as UTF-8 (BOM) CSVs so Chinese labels display correctly in Excel.

## Locale

Scripts use Chinese column labels in three-line tables. Run R under a UTF-8
locale (R ≥ 4.2 on Windows defaults to UTF-8); a `C` locale will fail to parse
Chinese identifiers.

## License

MIT — see [LICENSE](LICENSE).