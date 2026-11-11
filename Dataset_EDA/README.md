# Basic rclif_2_1 availability and overview EDA

Run from CLIF-Cohort on Discovery:

```bash
Rscript Dataset_EDA/eda_rclif_dataset.R ../rclif_2_1
```

Output defaults to a fresh `Dataset_EDA/results/basic_run_<timestamp>/` directory.
An optional second argument specifies another NEW output directory. Existing
runs are never overwritten, so older detailed results remain unchanged.

The updated script is intentionally basic. It profiles the delivered extract
before cohort construction. It does not read numerical lab/vital results or
medication doses, and creates no variable-by-variable distribution plots.

## Overview figures

At most four PNGs are created in plots/:

- sample_size_by_site.png: admissions and unique patients per hospital.
- admissions_over_time.png: monthly admissions overall and by hospital.
- table_sizes.png: number of rows in each source table (log scale).
- table_availability.png: percentage of source admissions with at least one row
  in each encounter-linked table.

Small site/month counts are omitted using SUMMARY_SMALL_CELL_THRESHOLD (default
10). A plot is skipped if it has no reportable data. Exact overall sample sizes,
table counts and record-availability denominators remain in the summaries.

## Basic CSV summaries

- sample_size_summary.csv: overall unique patients/admissions, duplicate and
  missing identifier counts, and conflicting admission-patient links.
- sample_size_by_site.csv: site-level admissions and unique patients.
- admissions_by_month.csv and admissions_by_site_and_month.csv.
- time_span.csv: admission/discharge coverage coarsened to first/last month,
  with valid and missing/invalid date counts; no exact timestamps exported.
- table_inventory.csv: source table names, row counts and number of fields.
- variable_dictionary.csv: every field name and Arrow data type.
- table_availability.csv: direct and linked unique-patient/admission counts,
  record coverage and orphan encounter-ID counts.
- measurement_types_available.csv: names of available lab, vital and numerical
  assessment categories only; no measurement values, counts or distributions.
- table_errors.csv: any incomplete inventory/availability steps.
- START_HERE.md and run_notes.txt: interpretation and run settings.

## Interpretation

An admission is a distinct hospitalization_id, not a merged hospital block.
Duplicate hospitalization rows are counted and deduplicated using earliest
recorded admission for the site/month summaries. Conflicting patient links are
flagged and excluded from linked-patient denominators. Overall direct patient
counts in the hospitalization table still reflect all distinct recorded IDs.

Hospital attribution uses the earliest valid recorded ADT hospital, with an
available hospitalization site field as fallback; otherwise the site is Unknown.
This is encounter attribution, not a claim of continuous ICU or hospital occupancy.

A patient may have multiple admissions and appear at multiple sites/months.
Therefore site/month patient counts should NOT be summed to infer overall unique
patients. The patient table may contain patients absent from the hospitalization
roster; these two sample sizes are reported separately.

Table availability means at least one row exists for an admission; it does not
establish complete measurements, low missingness within fields, or complete death
ascertainment. This simplified version does not profile individual-field missingness
or clinical-value distributions. Variable lists describe the supplied source, not
a validated clinical phenotype. Medication coverage is table-level only.

The extract's upstream ED/respiratory-support filters remain. Observed month
ranges do not establish continuous coverage or a data refresh date. Non-CLIF
shipping CSV/PDF files are not scanned.

## Runtime and output handling

Packages: arrow, data.table, dplyr, ggplot2. Tables are inventoried using schemas
and row counts; only identifier/category columns are batch-scanned for availability.
The hospitalization roster and exact identifier sets occupy memory proportional
to unique patients/admissions. Large vitals/medication value tables are not loaded.
Arrow and data.table run single-threaded.

Optional controls:

```bash
export EDA_BATCH_ROWS=65536
export EDA_MAX_UNIQUE_IDS=2000000
export SUMMARY_SMALL_CELL_THRESHOLD=10
```

The first limit also protects hospitalization-roster row size. Older numeric
reservoir/category-profile environment settings no longer apply. COMPLETE.txt
means all steps finished; PARTIAL.txt and table_errors.csv indicate missing
summaries. A missing/invalid core hospitalization roster stops execution.

No patient/encounter identifier VALUES, row-level observations, free-text values,
exact dates or detailed medication/measurement distributions are exported.
Identifiers remain transiently in memory for exact counting/linkage. These are
simpler working outputs, not automatic clearance to download or share; follow
normal project disclosure review.

Test with synthetic inputs only:

```bash
Rscript Dataset_EDA/tests/test_dataset_eda.R
```
