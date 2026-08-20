# mitox [![DOI](https://zenodo.org/badge/190498356.svg)](https://zenodo.org/badge/latestdoi/190498356)

Supportive scripts for the clinical trial, "A Pilot Study of the Effect of the Microbiome on Immune Checkpoint Inhibitor Response in Melanoma ([NCT05102773](https://clinicaltrials.gov/study/NCT05102773?term=spakowicz%20melanoma&rank=1), "OSU-19125)". The goal of the study is to predict which melanoma patients will respond to immunotherapy or experience immune-related adverse events using the gut microbiome and other biomarkers.

These analyses are reported in:
> Dravillas CE, Coleman SS IV, Hoyd R, Caryotakis G, Denko L, Chan CHF, et al. The Tumor Microbiome as a Predictor of Outcomes in Patients with Metastatic Melanoma Treated with Immune Checkpoint Inhibitors. Cancer Research Communications 2024;4:1978–90. https://doi.org/10.1158/2767-9764.CRC-23-0170.

Code to support the protocol submission, including power calculations, sample size estimates, estimated accrual dates, etc., can be found at:

> Spakowicz, D., and Muniak, M. (2019). spakowiczlab/mitox-protocol: IRB submission (Zenodo).

This work is supported by a Pelotonia (DS) and the Ohio State University's Translational Therapeutics Program (CV).

## Manuscript directory

Analyses for the manuscript live in `manuscript/`. Figure notebooks read derived objects from `manuscript/data/` and write images and table exports to `manuscript/figures/`.

| Path | Contents |
| --- | --- |
| `manuscript/scripts/processing/` | HPC batch scripts and the notebook that builds de-identified analysis objects |
| `manuscript/scripts/analysis/` | Statistical analyses and the notebooks that assemble each manuscript figure and table |
| `manuscript/data/` | Small lookup tables and model-comparison metrics committed with the repo |
| `manuscript/figures/` | Saved figure panels (SVG, and PNG when both were exported) plus table CSV/HTML exports |

Run notebooks from `manuscript/scripts/analysis/` (or `processing/`) so relative paths of the form `../../data/` and `../../figures/` resolve.

Clinical and microbiome objects used by the figure notebooks (for example `clinical.rds`, `counts-long.rds`, ANCOM-BC2 and random-forest outputs) are produced by the processing and analysis scripts below. They are not all stored in git.

### Processing scripts

| File | Role |
| --- | --- |
| `publication_objects.Rmd` | Builds the de-identified R objects that later notebooks load |
| `generate-metaphlan-batch.Rmd` | Writes the MetaPhlAn batch job used to profile shotgun reads |
| `merge_metaphlan-viruses.Rmd` | Aggregates MetaPhlAn taxon and virus tables |
| `humann-batch.Rmd` | Writes the HUMAnN functional-profiling batch job |

### Analysis scripts that feed the figures

These notebooks compute the results that the figure scripts plot. They do not write the manuscript figure files themselves.

| File | Role |
| --- | --- |
| `ANCOM-BC2.Rmd` | Baseline differential abundance for response and irAE |
| `ANCOM-BC2-respondersig.Rmd` | ANCOM-BC2 on the published Respondersig cohorts |
| `ANCOM-BC2-results-overlap.Rmd` | Overlap of response and irAE taxa |
| `respondersig.Rmd`, `respondersig-irae.Rmd` | Respondersig scores in this cohort |
| `alphadiversity.Rmd`, `alpha-diversity.Rmd` | Alpha-diversity summaries (Figure 2B / 3B) |
| `longitudinal-modeling-response.Rmd`, `longitudinal-modeling-irAE.Rmd` | Zero-inflated Gaussian mixed models (Figures 4 and 5) |
| `supplemental-figureS1.Rmd` | Longitudinal model comparison using `data/model-metrics.csv`; Supplemental Figure S1 |
| `core-microbiome.Rmd` | Per-patient core microbiome score (Figure 5E) |
| `random-forest.Rmd`, `rf-functions.R`, `corescore-predictive-modeling.Rmd`, `random-forest-validation.Rmd` | Random-forest models and external validation (Figure 6 and Supplemental Figures S3–S7) |
| `examine-timing.Rmd`, `window_size.Rmd` | Definitions of baseline, irAE-adjacent, and end-of-study samples |
| `labels_ANCOM-resist.Rmd`, `examine-viruses.Rmd` | Taxonomy labels and virus-table checks |

### Figures and the files that create them

Filenames in `manuscript/figures/` use the manuscript panel letter when that panel is a single exported plot. Multi-part schema panels and a few Illustrator composites are noted below.

#### Figure 1 — study schema and sampling timeline (`figure1.Rmd`)

| Panel | Saved file | Description |
| --- | --- | --- |
| 1A | `Fig1A-Baseline.svg` | Timing of baseline stool samples (within 10 days of C1D1) |
| 1A | `Fig1A-irAE.svg` | Timing of samples collected within 10 days of an irAE |
| 1A | `Fig1A-eos.svg` | Timing of ~12-week / end-of-study samples |
| 1C | `timeline-MET-restricted-dates.svg` | Per-patient event timeline for metastatic participants |

Figure 1B (CONSORT-style flow) is assembled outside these scripts.

#### Figure 2 — baseline microbiome and response (`figure2.Rmd`)

| Panel | Saved file | Description |
| --- | --- | --- |
| 2A | `Fig2A-response-baseline-stackedbar.svg` | Phylum stacked bars, baseline samples, faceted by response |
| 2B | `Fig2B-boxplots_diversity-response.svg` | Alpha diversity by response |
| 2C | `Fig2C-response-baseline-volcano.svg` | ANCOM-BC2 volcano plot, responders vs non-responders |
| 2D | `Fig2D-respondersig-plot.svg` | Respondersig ordination (written by the notebook; add to `figures/` after knitting) |
| 2E | `Fig2E-response-s0-baseline-boxplots.svg` | Structural-zero taxa at baseline by response |
| 2F | `Fig2F-rsig-overlap-effect-size.svg` | Effect sizes of taxa shared with Respondersig |

`Fig2E-labels.svg` is a legend/label companion for 2E.

#### Figure 3 — baseline microbiome and irAE (`figure3.Rmd`)

| Panel | Saved file | Description |
| --- | --- | --- |
| 3A | `Fig3A-irAE-baseline-stackedbar.svg` | Phylum stacked bars, baseline samples, faceted by irAE |
| 3B | `Fig3B_boxplots_diversity-irAE.svg` | Alpha diversity by irAE |
| 3C | `Fig3C-irae-baseline-volcano.svg` | ANCOM-BC2 volcano plot, irAE vs no irAE |
| 3D | `Fig3D-respondersig-plot.svg` | Respondersig ordination colored by irAE |
| 3E | `Fig3E-irAE-s0-baseline-boxplots.svg` | Structural-zero taxa at baseline by irAE |
| 3F | `Fig3F-pairwise-dist_boxplots.svg` | Pairwise distances between irAE and no-irAE groups |

`Fig3E-labels.svg` is a legend/label companion for 3E.

#### Figure 4 — longitudinal models of response (`figure4.Rmd`)

| Panel | Saved file | Description |
| --- | --- | --- |
| 4A | `Fig4A-ZIGMM-nonquad-volcano.svg` | ZIGMM volcano for the linear response × time term |
| 4B | `Fig4B-ZIGMM-quad-volcano.svg` | ZIGMM volcano for the quadratic time term |
| 4C | `Fig4C-lineplot.svg` | Longitudinal abundance of *Phocaeicola plebeius* by response |
| 4D | `Fig4D-venn.svg` | Overlap of ANCOM-BC2 and ZIGMM taxa |

The notebook currently also writes PNG copies of 4A–4C.

#### Figure 5 — longitudinal models of irAE (`figure5.Rmd`)

| Panel | Saved file | Description |
| --- | --- | --- |
| 5A | `Fig5A-ZIGMM-nonquad-volcano-irae.svg` | ZIGMM volcano for irAE-associated taxa |
| 5B | `Fig5C-lineplot.svg` | Longitudinal abundance of *Phocaeicola dorei*, with irAE event times marked |
| 5C | `Fig5C-venn.svg` | Overlap of ANCOM-BC2 and ZIGMM irAE taxa |
| 5E | `Fig5E-boxplot.svg` | Core microbiome score by irAE |

Figure 5D (schema) is assembled outside these scripts.

#### Figure 6 — prediction of response and irAE (`figure6.Rmd`)

| Panel | Saved file | Description |
| --- | --- | --- |
| 6A | `Fig6A-Forestplot-response.svg` | Univariate clinical predictors of response |
| 6B | `Fig6B-Forestplot-irae.svg` | Univariate clinical predictors of irAE |
| 6C | `Fig6C-AUCPR-boxplots-response.svg` | Random-forest AUPRC by feature set, response |
| 6D | `Fig6D-AUCPR-boxplots-irae.svg` | Random-forest AUPRC by feature set, irAE |
| 6E | `Fig6E-corescore-AUCPR-boxplots-irae.svg` | AUPRC when the core score is included, irAE |
| — | `Fig6I-combined.svg` | External-validation ROC curves by published cohort |
| — | `Fig6J-varImp-response.svg` | Variable importance for response models |
| — | `Fig6K-varImp-irae.svg` | Variable importance for irAE models |

`Fig6I-valROC-by-study.svg` is an alternate validation plot (not used in the main figure). Earlier ROC exports remain in `figures/` as `Fig6E-ROC-response.svg`, `Fig6F-ROC-irae.svg`, `Fig6G-valROC-response.svg`, and `Fig6H-valROC-irae.svg`; those views now correspond to Supplemental Figures S3, S4, S6, and S7. `Fig6A-Forestplot-irae.svg` and `Fig6B-Forestplot-response.svg` are previous filename assignments (the current 6A/6B files above are the ones used in the manuscript). `corescore-AUCPR-boxplots-response.svg` and `corescore-AUCPR-boxplots-irae.svg` are exploratory core-score plots.

#### Supplemental figures

| Figure | Script | Saved file | Description |
| --- | --- | --- | --- |
| S1 | `supplemental-figureS1.Rmd` | `supplemental-model-metrics.svg` | Comparison of longitudinal model families (AIC/MSE); also writes `data/model-metrics.csv` |
| S2A | `supplemental-figureS2-S7.Rmd` | `supplemental-figuresS2A.svg` | Full univariate clinical forest plot, response |
| S2B | `supplemental-figureS2-S7.Rmd` | `supplemental-figureS2B.svg` | Full univariate clinical forest plot, irAE |
| S3 | `supplemental-figureS2-S7.Rmd` | `supplemental-figureS3.svg` | ROC: microbiome vs clinical models, response |
| S4 | `supplemental-figureS2-S7.Rmd` | `supplemental-figureS4.svg` | ROC: microbiome vs clinical models, irAE |
| S5 | `supplemental-figureS2-S7.Rmd` | `supplemental-figureS5.png` | Variable importance for irAE models that include the core score |
| S6 | `supplemental-figureS2-S7.Rmd` | `supplemental-figureS6.png` | External-validation ROC for response (all cohorts combined) |
| S7 | `supplemental-figureS2-S7.Rmd` | `supplemental-figureS7.png` | External-validation ROC for irAE |

S5–S7 are written as PNG only. `window_size.Rmd` writes diagnostic plots `irae_windowsize.png` and `irae_windowsize_closeup.png` that are not manuscript figures.

### Tables

| Table | Script | Saved file | Description |
| --- | --- | --- | --- |
| 1 | `tables.Rmd` | `table1.html` | Cohort characteristics stratified by response |
| 2 | `tables.Rmd` | `table2.html` | Cohort characteristics stratified by irAE |
| 3 | `tables.Rmd` | `table3.html` | irAE types and grades |
| S1 | `tables.Rmd` | `supplemental-table1.csv` | irAE type/grade counts (CSV export of Table 3) |
| S2A/B | `supplemental-tables.Rmd` | `supplemental-table2A.csv`, `supplemental-table2B.csv` | ANCOM-BC2 response results and structural zeros |
| S3A/B | `supplemental-tables.Rmd` | `supplemental-table3A.csv`, `supplemental-table3B.csv` | ANCOM-BC2 Respondersig results and structural zeros |
| S4A/B | `supplemental-tables.Rmd` | `supplemental-table4A.csv`, `supplemental-table4B.csv` | ANCOM-BC2 irAE results and structural zeros |
| S5 | `supplemental-tables.Rmd` | `supplemental-table5.csv` | ZIGMM coefficients for response |
| S6 | `supplemental-tables.Rmd` | `supplemental-table6.csv` | ZIGMM coefficients for irAE |
| S7 | `supplemental-tables.Rmd` | `supplemental-table7.csv` | Participant-level clinical variables |
| S8 | `supplemental-tables.Rmd` | `supplemental-table8.csv` | Stool-sample metadata |

### Data files in `manuscript/data/`

| File | Description |
| --- | --- |
| `model-metrics.csv` | Fit statistics comparing longitudinal model families; input to Supplemental Figure S1 |
| `CRISPR_VSG-to-species.csv` | Virus sequence group to species lookup |
| `CRISPR_VSG-to-SGBs.csv` | Virus sequence group to SGB lookup | 
