# adjing [![DOI](https://zenodo.org/badge/651603335.svg)](https://doi.org/10.5281/zenodo.21841704)
Scripts to regenerate all analyses and figures presented in:
> ## Pre-treatment T-cell Transcriptional Signatures Predict Immunotherapy Outcomes in Melanoma
> *Noah Lepola<sup>1,2</sup>, Caroline Dravillas<sup>3</sup>, Shannon Gray<sup>1,2</sup>, Michael S. Bodnar<sup>1,2</sup>, Namrata Arya<sup>1,2</sup>, Richard Wu<sup>3</sup>, Claire Verschraegen<sup>3</sup>, William E Carson<sup>3</sup>, Kari L Kendra<sup>3</sup>, Daniel J Spakowicz<sup>3,4^</sup>, Christin E Burd<sup>1,2,4^</sup>*<br>
><sup>1</sup> Department of Molecular Genetics, The Ohio State University College of Arts and Sciences, Columbus, Ohio<br>
><sup>2</sup> Department of Cancer Biology and Genetics, The Ohio State University College of Medicine, Columbus, Ohio<br>
><sup>3</sup> Division of Medical Oncology, Department of Internal Medicine, The Ohio State University Comprehensive Cancer Center; Columbus, OH, USA.<br>
><sup>4</sup> Pelotonia Institute for Immuno-Oncology, The Ohio State University Comprehensive Cancer Center; Columbus, OH, USA.<br>

submitted to *bioRxiv*

## Repository structure

| Path | Contents |
| --- | --- |
| `manuscript/` | Scripts, shared helpers, and outputs used to regenerate the manuscript figures and supporting tables |
| `manuscript/scripts/` | R Markdown notebooks that process data and create each figure panel |
| `manuscript/figures/` | PNG figure panels written by the scripts above |
| `manuscript/Outputs/` | Intermediate R objects (`.RData` / `.rds`), model fits, and derived tables consumed by the figure scripts |
| `manuscript/rf_functions.R` | Shared random-forest and plotting helpers sourced by the figure notebooks |
| `exploratory/` | Earlier exploratory analyses not used for the final manuscript figures |
| `00-paths.R` / `adjing.json` | Local path configuration for clinical and NanoString source data |
| `adjing.Rproj` | RStudio project file |

The `manuscript/` directory is the entry point for regenerating the published analyses. Numbered notebooks (`00_`, `01_`, `02_`) prepare adjuvant and metastatic cohort objects; subsequent notebooks train models and write figure panels into `manuscript/figures/`. Script filenames encode the figure panels they produce (for example, `adj_rec_Fig2ABCDGH.Rmd` generates Figure 2 panels A–D, G, and H).

### Figure-to-script map

Use this table to locate the notebook that regenerates a given panel. Output PNGs are written under `manuscript/figures/` unless noted otherwise.

#### Main figures

| Figure panels | Script | Output files |
| --- | --- | --- |
| Fig. 2A, 2B, 2C, 2D, 2G, 2H | [`adj_rec_Fig2ABCDGH.Rmd`](manuscript/scripts/adj_rec_Fig2ABCDGH.Rmd) | `AdjRecLogr.png` (2A), `AdjRecClinical.png` (2B), `AdjRecFullpanel.png` (2C), `AdjRecBoruta.png` (2D), `AdjRecDistribution.png` (2H) |
| Fig. 2E | [`logr_adj_rec_and_volcano_Fig_2E.Rmd`](manuscript/scripts/logr_adj_rec_and_volcano_Fig_2E.Rmd) | `RecAdjVol.png` |
| Fig. 2F, Fig. 3D | [`F1_opt_adj_Fig2F_Fig3D.Rmd`](manuscript/scripts/F1_opt_adj_Fig2F_Fig3D.Rmd) | `RecAdjOptBar.png` (2F), `ToxAdjOptBar.png` (3D) |
| Fig. 3A, 3B, 3E | [`adj_tox_Fig3ABE.Rmd`](manuscript/scripts/adj_tox_Fig3ABE.Rmd) | `AdjToxDualGraph.png` (3A), `AdjToxBoruta.png` (3B), `AdjToxLogr.png` (3E) |
| Fig. 3C | [`logr_adj_tox_and_volcano_Fig_3C.Rmd`](manuscript/scripts/logr_adj_tox_and_volcano_Fig_3C.Rmd) | `ToxAdjVol.png` |
| Fig. 5A, 5B, 5E | [`meta_pro_Fig5ABE.Rmd`](manuscript/scripts/meta_pro_Fig5ABE.Rmd) | `MetaProDualGraph.png` (5A), `MetaProBoruta.png` (5B), `MetaProLogr.png` (5E) |
| Fig. 5C | [`logr_meta_pro_and_volcano_Fig_5C.Rmd`](manuscript/scripts/logr_meta_pro_and_volcano_Fig_5C.Rmd) | `ProMetVol.png` |
| Fig. 5D, Fig. S1D | [`F1_opt_meta_Fig5D_FigS1D.Rmd`](manuscript/scripts/F1_opt_meta_Fig5D_FigS1D.Rmd) | `ProMetOptBar.png` (5D), `ToxMetOptBar.png` (S1D) |
| Fig. 6A, 6B; Fig. S2A, S2B | [`cross_setting_comparison_Fig6AB_FigS2AB.Rmd`](manuscript/scripts/cross_setting_comparison_Fig6AB_FigS2AB.Rmd) | `AdjTrainRecCrossCohort.png` (6A), `MetaTrainRecCrossCohort.png` (6B), `AdjTrainToxCrossCohort.png` (S2A), `MetaTrainToxCrossCohort.png` (S2B) |

#### Supplementary figures

| Figure panels | Script | Output files |
| --- | --- | --- |
| Fig. S1A, S1B, S1E | [`meta_tox_FigS1ABE.Rmd`](manuscript/scripts/meta_tox_FigS1ABE.Rmd) | `MetaToxDualGraph.png` (S1A), `MetaToxBoruta.png` (S1B), `MetaToxLogr.png` (S1E) |
| Fig. S1C | [`logr_meta_tox_and_volcano_Fig_S1C.Rmd`](manuscript/scripts/logr_meta_tox_and_volcano_Fig_S1C.Rmd) | `ToxMetVol.png` |
| Fig. S1D | [`F1_opt_meta_Fig5D_FigS1D.Rmd`](manuscript/scripts/F1_opt_meta_Fig5D_FigS1D.Rmd) | `ToxMetOptBar.png` |
| Fig. S2A, S2B | [`cross_setting_comparison_Fig6AB_FigS2AB.Rmd`](manuscript/scripts/cross_setting_comparison_Fig6AB_FigS2AB.Rmd) | `AdjTrainToxCrossCohort.png` (S2A), `MetaTrainToxCrossCohort.png` (S2B) |

#### Data preparation and supporting tables

These notebooks do not write main-figure panels, but they create the cohort objects, clinical tables, and intermediate results used by the figure scripts:

| Script | Role |
| --- | --- |
| [`00_transcript_processing_adjuvant.Rmd`](manuscript/scripts/00_transcript_processing_adjuvant.Rmd) | Process adjuvant NanoString / clinical characteristics |
| [`01_outcome_processing_adjuvant.Rmd`](manuscript/scripts/01_outcome_processing_adjuvant.Rmd) | Process adjuvant patient outcomes |
| [`02_object_processing_metastatic.Rmd`](manuscript/scripts/02_object_processing_metastatic.Rmd) | Build metastatic cohort clinical and outcome objects |
| [`tableone_adjuvant.Rmd`](manuscript/scripts/tableone_adjuvant.Rmd) | Adjuvant cohort Table 1 |
| [`tableone_metatatic.Rmd`](manuscript/scripts/tableone_metatatic.Rmd) | Metastatic cohort Table 1 |
| [`boruta_transcript_counts.Rmd`](manuscript/scripts/boruta_transcript_counts.Rmd) | Boruta feature-selection transcript counts across outcomes |
| [`adjing timeline.Rmd`](manuscript/scripts/adjing timeline.Rmd) | Adjuvant cohort timelines (`manuscript/Outputs/`) |
| [`metastatic timeline.Rmd`](manuscript/scripts/metastatic timeline.Rmd) | Metastatic cohort timelines (`manuscript/Outputs/`) |

To regenerate a figure, open the corresponding notebook under `manuscript/scripts/`, ensure paths in `adjing.json` point to the source data on your system, knit the upstream processing notebooks if `manuscript/Outputs/` objects are not already present, then knit the figure notebook.
