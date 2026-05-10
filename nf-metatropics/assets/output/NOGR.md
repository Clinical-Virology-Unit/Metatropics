# NOGR (Non-Overlapping Genomic Regions)

Metatropics runs viral classification with Virasign. NOGR (non-overlapping genomic regions) is computed in that Virasign step for each candidate virus from reads aligned to the selected reference, then carried into the Metatropics summary report (HTML and CSV) next to breadth, depth, and related columns. It summarizes how many **distinct, non-overlapping stretches** of the reference genome have support from at least one mapped read, and how many reference bases those stretches cover in total. That matches a recurring mNGS idea: insist on several **separated** regions of evidence so a narrow pile-up (for example amplicon contamination) is less likely to pass as a genome-wide signal—see Nature Communications [`s41467-024-51470-y`](https://www.nature.com/articles/s41467-024-51470-y) for one published rule (≥3 non-overlapping viral reads/contigs as a positivity aid).

---

## Why this is useful

Coverage breadth can be low for real detections (degraded libraries, low viral load) and also for artifacts where many reads stack on one small region of the genome. NOGR adds a second axis: how spread out the mapped evidence is along the reference. That helps separate plausible genome-spanning signal from tight technical noise. Conversely, when coverage is already very high and near-complete, mapped intervals merge and the non-overlapping region count is usually small, so breadth and depth already argue strongly for the virus and NOGR adds little beyond that.

---

## What exactly is computed

Virasign uses the per-virus BAM from classification (reads aligned to that virus’s chosen reference). Each mapped read contributes its interval on the reference; overlapping intervals are merged so you end up with a set of disjoint regions that still have read support.

Two numbers are stored with each hit:

- `nogr_regions`: how many disjoint regions that is
- `nogr_bases`: how many reference bases lie inside those regions together (never more than the reference length)

In the Metatropics summary HTML and CSV, both numbers appear in one column headed **NOGR (#/bases)**, formatted as `regions|bases`.

---

## Where it appears

You can read NOGR in `Summary/metatropics/`, for example `Metatropics_Summary_RVDB.html` and `Metatropics_Summary_RVDB.csv`, in the NOGR (#/bases) column—the same heading in HTML and CSV. A minimum region count filter, when you want it, is configured with `--virasign_min_nogr`; see [`../submission/all_options.md`](../submission/all_options.md).

For Virasign-only details, see [`DaanJansen94/virasign`](https://github.com/DaanJansen94/virasign).

---
