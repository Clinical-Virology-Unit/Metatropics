# Z-score background correction

Metatropics runs viral classification with Virasign. A background-corrected Z-score per virus is computed in that Virasign step from your negative controls (for example water or extraction blanks), then incorporated into the Metatropics summary report (HTML and CSV) next to each detection. This mirrors a common idea in metagenomic reporting: quantify whether a signal is unusually high compared to background contamination.

---

## Why use a Z-score?

Some taxa show up at low levels in many runs due to:
- reagent / extraction contamination
- index hopping / low-level carryover
- laboratory environment background

By comparing each virus to your negative controls in the same run, the Z-score answers:

> “Is this virus higher than what we typically see in negative controls?”

---

## When it is computed in a Metatropics run

In a Metatropics run, Z-scores are only produced when at least two negative controls are present; with fewer than two, no Z-score is computed. That matches how the score is defined: it expresses how far a sample lies from the mean of the controls in units of the controls’ standard deviation (see [Formula](#formula)). A single control fixes a reference level but does not define a spread across the background, so the “standard deviation units” part of the Z-score is not meaningful until there are at least two negative controls.

When the minimum of two negative controls is met, auto-detection is the default: those controls are inferred from sample names that contain `water`, `h2o`, or `h20` (case-insensitive), reflecting common naming for blanks. How to configure this in a run is listed with the other pipeline options in [`../submission/all_options.md`](../submission/all_options.md).

---

## What signal is used

The Z-score is computed from the mapped read count for each virus: how many reads from the sample align to that virus’s reference. That is the same mapped-read count shown for that detection in the Metatropics summary report.

---

## Formula

For each virus, the Z-score is the number of standard deviations that the sample’s log-transformed mapped_reads lies above or below the mean of the log-transformed mapped_reads in the selected negative controls.

This follows the same background-correction idea used by CZ ID / IDseq background models. For more information, see:

- IDseq paper (GigaScience, 2020): `https://ncbi.nlm.nih.gov/pmc/articles/PMC7566497/`
- CZ ID workflows wiki: `https://github.com/chanzuckerberg/czid-workflows/wiki`

---

## How to interpret Z-scores

Z-scores are “standard deviation units” above/below background. Practical interpretation depends on context and other evidence (coverage breadth/depth, identity, NOGR, etc.), but the table below is a useful starting point.

| Z-score | SD interpretation | Practical interpretation |
|---:|---|---|
| -3 | 3 SD below mean | Below background |
| -1 | 1 SD below mean | Background-like |
| 0 | At mean | Background-like |
| 1 | 1 SD above mean | Slightly above background |
| 2 | 2 SD above mean | Above background |
| 3 | 3 SD above mean | Strongly above background |
| 10 | 10 SD above mean | Strongly above background |
| 50 | Extreme outlier | Strongly above background |
| 100 | Extreme outlier | Only in samples that are not negative controls

---

## Where it appears in Metatropics outputs

Under `Summary/metatropics/` (for example `Metatropics_Summary_RVDB.html` and `Metatropics_Summary_RVDB.csv`), the metric appears as a column titled **Z-score** in both the HTML report and the CSV export—not only in one of them. Each detection row carries a value in that column. Only in the CSV do you also get the **Sample** and **Background** columns: Background is `yes` on negative-control samples used for the Z-score background.

For Virasign-only details or running the classifier outside Metatropics, see [`DaanJansen94/virasign`](https://github.com/DaanJansen94/virasign).
