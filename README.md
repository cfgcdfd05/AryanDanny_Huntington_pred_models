# Huntington's Disease m6A Biomarker Classification

This project investigates whether a panel of 12 m6A (RNA methylation) related genes can distinguish Huntington's disease (HD) samples from non-HD samples, using gene expression data pulled from four independent public datasets on the Gene Expression Omnibus (GEO). Each dataset is extracted into a common ML-ready CSV format and used to train the same set of four classifiers (SVM, Logistic Regression, Random Forest, XGBoost).

## Biomarker panel

All four datasets are restricted to the same 12 candidate genes:

```
TUT7, PPP1CC, CEBPB, CEBPD, FTO, METTL16,
IGF2BP3, YTHDF1, YTHDF2, YTHDF3, YTHDC1, YTHDC2
```

These genes are all involved in writing, reading, or erasing the m6A RNA modification. Not every gene is present in every dataset — coverage depends on the microarray platform used or on upstream filtering during RNA-seq quantification (see per-dataset notes below).

---

## Data sources

| Dataset | Platform | Tissue | Groups | GEO link |
|---|---|---|---|---|
| GSE3790 | Affymetrix HG-U133A/B (GPL96 + GPL97) | — | HD vs. control | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3790 |
| GSE1751 | Affymetrix HG-U133A (GPL96) | Whole blood (PAXgene) | Symptomatic HD, presymptomatic HD, healthy control | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1751 |
| GSE64810 | Illumina HiSeq 2000 (RNA-seq) | Brain, BA9 prefrontal cortex | HD vs. neurologically normal control | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810 |
| GSE33000 | Custom Agilent array (GPL4372) | Brain, prefrontal cortex | Huntington's disease, Alzheimer's disease, non-demented control | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33000 |

### GSE3790
Two-platform Affymetrix microarray series (GPL96 + GPL97). Sample labels are parsed from the `!Sample_characteristics_ch1` / `!Sample_title` lines of each platform's series matrix file, using keyword matching for `hd`/`huntington` (label 1) vs. `control`/`normal` (label 0). Probe-to-gene mapping comes from each platform's `.annot` file. Data from both platforms is merged per `GSM` sample ID before probe averaging.

### GSE1751
Single-platform Affymetrix series (GPL96), whole blood samples from a Krainc lab study (Borovecki et al.). Labels come from the `!Sample_description` line, which contains free text like `"Huntington's disease patient 1, symptomatic"` or `"healthy control 1"`. This dataset also contains 5 **presymptomatic** HD gene-carriers (`P1`-`P5`) — confirmed HD mutation carriers with no current symptoms. These are excluded from the final dataset to keep a clean binary symptomatic-vs-control comparison (12 HD, 14 control, 26 total samples).

### GSE64810
RNA-seq dataset (BA9 prefrontal cortex), quantified as DESeq2-normalized counts, indexed by Ensembl gene ID (e.g. `ENSG00000083223.10`) rather than probe ID. Sample IDs follow an `H_xxxx` (HD) / `C_xxxx` (control) naming convention, which is used directly for labeling. Because this is RNA-seq, each gene is already a single row — no probe averaging is needed. Two of the 12 target genes, **CEBPD** and **YTHDF2**, are absent from this particular count matrix (likely filtered out upstream during DESeq2 processing due to low expression), so this dataset trains on 10 of 12 genes.

### GSE33000
Custom Agilent array (GPL4372), prefrontal cortex tissue, containing three disease-status groups: Alzheimer's disease, Huntington's disease, and non-demented control (parsed from a `!Sample_characteristics_ch2` line starting with `"disease status:"`). Alzheimer's disease samples are excluded during extraction to preserve a binary HD-vs-control setup consistent with the other three datasets.

---

## Dataset CSV format

All four output CSVs share the same structure:

| Column | Description |
|---|---|
| `Sample_ID` | GEO sample accession (`GSM...`) or cohort-specific sample name (e.g. `C_0002`, `H_0001`) |
| `Huntingtons_Disease` | Target label: `1` = Huntington's disease, `0` = control |
| one column per gene | Expression value for that gene, averaged across all probes mapping to it (microarray datasets) or taken directly from the normalized count matrix (RNA-seq dataset) |

Rows are **samples** (patients), columns are **gene features** plus the target — the standard orientation for scikit-learn (`X` = feature matrix, `y` = target vector), obtained by transposing the raw GEO series matrix (which stores genes as rows and samples as columns).

For microarray datasets, a single gene is often measured by multiple probes on the same array (e.g. `IGF2BP3_203819_s_at` and `IGF2BP3_203820_s_at`). During extraction, these are temporarily kept as separate `GeneSymbol_ProbeID` columns, then averaged row-wise into one column per gene symbol in the final CSV, since a model trained directly on the CSV would otherwise be unable to distinguish probes intended to measure the same underlying transcript.

---

## Models

Each notebook trains and compares four classifiers, using `GridSearchCV` to tune hyperparameters and `StratifiedKFold` cross-validation (fold count automatically capped by the minority class size, since these datasets are small) to select the best configuration for each.

### SVM (Support Vector Machine)
Finds the decision boundary that maximizes the margin between HD and control samples in feature space. Tuned over kernel type (`linear` vs. `rbf`) and the regularization strength `C`. Probability outputs (needed for ROC/AUC) come from an internal calibration step rather than the raw decision function.

### GLM (Logistic Regression)
A linear model that estimates the probability of HD status as a sigmoid-transformed weighted sum of gene expression values. Tuned over the regularization strength `C`. Coefficients are directly interpretable (a positive coefficient means higher expression pushes the prediction toward HD).

### Random Forest
An ensemble of decision trees, each trained on a bootstrap sample of the data using a random subset of features per split, with predictions averaged across trees. Tuned over number of trees, max depth, and minimum samples required to split/form a leaf. Captures non-linear relationships and feature interactions that logistic regression cannot.

### XGBoost
A gradient-boosted ensemble of decision trees, where each new tree is fit to correct the errors of the previous trees (as opposed to Random Forest's independent trees). Tuned over number of trees, max depth, learning rate, and subsample ratio. Generally the strongest performer on structured/tabular data of this kind, though also the most prone to overfitting on very small sample sizes.

### Evaluation
For each model: ROC curve and AUC (on a held-out test split), confusion matrix, accuracy, and a full classification report (precision/recall/F1 per class). Random Forest and XGBoost additionally report feature importances, showing which of the 12 genes contributed most to each model's predictions.

Given the small sample sizes in these datasets (ranging from roughly 13 to a few hundred samples depending on the cohort), AUC values from a single train/test split should be read cautiously — with few test samples, AUC can only take a limited set of discrete values and is sensitive to how individual samples happen to fall into the test set. Cross-validated evaluation (as used inside `GridSearchCV`, and in the per-gene univariate analysis in later notebook sections) gives a more stable picture than any single test-split number.

---

## Caveats and known limitations

- **Small sample sizes.** All four datasets are small by machine learning standards (tens to low hundreds of samples). Results should be interpreted as exploratory rather than clinically validated.
- **GSE64810 is missing 2 of 12 genes** (CEBPD, YTHDF2) due to upstream filtering in the provided normalized count matrix — models for this dataset are trained on 10 genes, not 12.
- **GSE1751 excludes 5 presymptomatic HD samples** to preserve a clean binary comparison; these carry the HD mutation but have not yet developed symptoms.
- **GSE33000 excludes all Alzheimer's disease samples**, keeping only HD vs. non-demented control.
- **Tissue and platform differ across datasets** (blood vs. brain; microarray vs. RNA-seq), so a gene's importance in one dataset does not necessarily generalize to another — genes that rank highly for HD prediction consistently across multiple datasets are a stronger signal than agreement within just one.
